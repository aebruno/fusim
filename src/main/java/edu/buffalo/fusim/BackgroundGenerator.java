/*
 * Copyright 2012 Andrew E. Bruno <aebruno2@buffalo.edu>
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not
 * use this file except in compliance with the License. You may obtain a copy
 * of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
 * License for the specific language governing permissions and limitations
 * under the License.
 */

package edu.buffalo.fusim;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.TimeUnit;

import net.sf.samtools.AbstractBAMFileIndex;
import net.sf.samtools.BAMIndexMetaData;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import cern.colt.list.IntArrayList;
import cern.jet.random.sampling.RandomSamplingAssistant;
import edu.buffalo.fusim.gtf.GTFParseException;

public class BackgroundGenerator implements FusionGenerator {
    private Log logger = LogFactory.getLog(BackgroundGenerator.class);

    private GeneModelParser parser;
    private ArrayBlockingQueue<TranscriptRecord> queue;
    private File backgroundFile;
    private double rpkmCutoff;
    private int threads;
    private boolean useBins;

    public BackgroundGenerator(File backgroundFile, GeneModelParser parser, double rpkmCutoff, int threads, boolean useBins) {
        this.parser = parser;
        this.queue = new ArrayBlockingQueue<TranscriptRecord>(100000);
        this.backgroundFile = backgroundFile;
        this.rpkmCutoff = rpkmCutoff;
        this.threads = threads;
        this.useBins = useBins;
    }

    public List<FusionGene> generate(File gtfFile, int nFusions) {
        logger.info("Processing background reads...");
        logger.info("Computing RPKM values using " + threads + " threads...");

        GeneModelProducer producer = new GeneModelProducer(gtfFile);
        producer.start();

        // Leave one thread for the GTF producer
        if (threads > 1)
            threads--;
        ArrayList<GeneModelConsumer> consumers = new ArrayList<GeneModelConsumer>();

        for (int i = 0; i < threads; i++) {
            GeneModelConsumer consumer = new GeneModelConsumer(backgroundFile, this.rpkmCutoff);
            consumer.start();
            consumers.add(consumer);
        }

        try {
            producer.join();
            for (GeneModelConsumer c : consumers) {
                c.join();
            }
        } catch (InterruptedException e) {}
        
        List<RPKM> rpkm = new ArrayList<RPKM>();
        
        // Reduce
        for(GeneModelConsumer c : consumers) {
            for(RPKM r : c.getRPKMList()) {
                rpkm.add(r);
            }
        }
        
        Collections.sort(rpkm, new RPKMCompare());
        
        List<FusionGene> fusions = new ArrayList<FusionGene>();

        // First bin genes into RPKM buckets
        BinGenes bin = new BinGenes(nFusions);
        bin.fill(rpkm);

        if(useBins) {
            logger.info("Generating fusions using binned RPKM values...");
            for(int i = 0; i < bin.size(); i++) {
                IntArrayList b = bin.getBin(i);
                b.trimToSize();
                int[] sample = RandomSamplingAssistant.sampleArray(2, b.elements());
                fusions.add(new FusionGene(rpkm.get(sample[0]).getTranscript(),
                                           rpkm.get(sample[1]).getTranscript()));
            }
        } else {
            logger.info("Generating fusions based on uniform background distribution...");
            Random r = new Random();
            int[] distribution = new int[rpkm.size()];
            int index = 0;
            for(int i = 0; i < bin.size(); i++) {
                IntArrayList b = bin.getBin(i);
                b.trimToSize();
                for(int j = 0; j < b.size(); j++) {
                    distribution[index] = i;    
                    index++;
                }
            }

            for(int x = 0; x < nFusions; x++) {
                int binIndex1 = distribution[r.nextInt(distribution.length)];
                int binIndex2 = distribution[r.nextInt(distribution.length)];
                IntArrayList b1 = bin.getBin(binIndex1);
                IntArrayList b2 = bin.getBin(binIndex2);

                fusions.add(new FusionGene(rpkm.get(b1.get(r.nextInt(b1.size()))).getTranscript(),
                                           rpkm.get(b2.get(r.nextInt(b2.size()))).getTranscript()));
            }
        }
        
        return fusions;
    }
  
    protected class RPKM {
        private TranscriptRecord transcript;
        private double value;
        
        public RPKM(TranscriptRecord transcript, double value) {
            this.transcript = transcript;
            this.value = value;
        }

        public TranscriptRecord getTranscript() {
            return transcript;
        }

        public double getValue() {
            return value;
        }
    }
    
    protected class RPKMCompare implements Comparator<RPKM> {
        public int compare(RPKM o1, RPKM o2) {
            return Double.compare(o1.getValue(), o2.getValue());
        }
    }
    
    protected class BinGenes {
        private IntArrayList[] bins;
        private int binsize;
        
        public BinGenes(int binsize) {
            this.bins = new IntArrayList[binsize];
            for(int i = 0; i < binsize; i++) {
                bins[i] = new IntArrayList();
            }
            this.binsize = binsize;
        }
        
        public void fill(List<RPKM> values) {
            //XXX fix this.. wtf. note values needs to be sorted
            int inc = (int)Math.ceil(values.size()/binsize);
            int binI = 0;
            for(int i = 0; i < values.size(); i++) {
                if((i+1)%inc == 0) {
                    binI++;
                    if(binI > binsize-1) binI = binsize-1;
                }
                bins[binI].add(i);
            }
        }
        
        public int size() {
            return this.binsize;
        }
        
        public IntArrayList getBin(int index) {
            return this.bins[index];
        }
        
    }

    protected class GeneModelProducer extends Thread {
        private Log log = LogFactory.getLog(GeneModelProducer.class);

        private File gtfFile;

        public GeneModelProducer(File gtfFile) {
            this.gtfFile = gtfFile;
        }

        public void run() {
            try {
                Charset charset = Charset.forName("UTF-8");

                // XXX Use Java NIO for reading File. Requires Java 7
                BufferedReader reader = Files.newBufferedReader(FileSystems
                        .getDefault().getPath(gtfFile.getAbsolutePath()),
                        charset);

                String line = null;
                while ((line = reader.readLine()) != null) {
                    if (line.startsWith("#")) continue;
                    TranscriptRecord feature;
                    try {
                        feature = parser.parseLine(line);
                        if(feature == null) continue;

                        //XXX skip the haplotypes and unassembled chroms
                        if(feature.getChrom().contains("_")) continue;
                        queue.put(feature);
                    } catch (InterruptedException e) {
                        log.fatal("InterruptedException while adding to queue: "
                                + e.getMessage());
                    } catch(GTFParseException e) {
                        log.fatal("Failed to parse gene model line: "
                                + e.getMessage());
                    }
                }
            } catch (IOException e) {
                log.fatal("I/O error while reading GTF file: "
                        + e.getMessage());
            }
        }
    }

    protected class GeneModelConsumer extends Thread {
        private Log log = LogFactory.getLog(GeneModelProducer.class);
        private SAMFileReader sam;
        private int totalMappedReads;
        private double rpkmCutoff;
        private List<RPKM> rpkmList = new ArrayList<RPKM>();

        public GeneModelConsumer(File bamFile, double rpkmCutoff) {
            this.rpkmCutoff = rpkmCutoff;
            File bamIndexFile = new File(bamFile.getAbsolutePath() + ".bai");
            if (bamIndexFile.canRead()) {
                sam = new SAMFileReader(bamFile, bamIndexFile);
            } else {
                //XXX this needs to throw an error!!! we need a BAM index
                sam = new SAMFileReader(bamFile);
            }

            AbstractBAMFileIndex index = (AbstractBAMFileIndex) sam.getIndex();

            totalMappedReads = 0;
            for (int i = 0; i < index.getNumberOfReferences(); i++) {
                BAMIndexMetaData meta = index.getMetaData(i);
                totalMappedReads += meta.getAlignedRecordCount();
            }

            //logger.warn("Total Mapped Reads found: "+totalMappedReads);
        }
        
        public List<RPKM> getRPKMList() {
            return this.rpkmList;
        }

        public void run() {
            TranscriptRecord transcript;

            while (true) {
                try {
                    transcript = queue.poll(100, TimeUnit.MILLISECONDS);
                } catch (InterruptedException e) {
                    log.fatal("Interrupted exception while consuming: "
                            + e.getMessage());
                    return;
                }
                if (transcript == null)
                    break;

                int count = 0;
                // XXX do we want only coding exons here???
                //for(int[] exon : feature.getCodingExons()) 
                for(int i = 0; i < transcript.getExonStarts().length; i++) {
                    // XXX end-1 here???
                    int start = transcript.getExonStarts()[i];
                    int end = transcript.getExonEnds()[i];
                    SAMRecordIterator it = sam.queryOverlapping(transcript.getChrom(), start, end);

                    while (it.hasNext()) {
                        SAMRecord samRecord = it.next();
                        // XXX do we require the mate to be mapped??
                        //if(samRecord.getReadUnmappedFlag() || samRecord.getMateUnmappedFlag()) continue;
                        //if(samRecord.getReadUnmappedFlag()) continue;
                        if(!samRecord.getReadUnmappedFlag() || !samRecord.getMateUnmappedFlag()) {
                            count++;
                        }
                    }

                    it.close();
                }
                Double rpkm = (Math.pow(10,9)*(double)count)/((double)totalMappedReads*transcript.getExonBases());
                if(rpkm > this.rpkmCutoff) {
                    rpkmList.add(new RPKM(transcript, rpkm));
                }
            }
        }
    }
}
