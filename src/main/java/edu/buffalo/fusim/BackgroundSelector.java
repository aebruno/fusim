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
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.TimeUnit;

import net.sf.samtools.AbstractBAMFileIndex;
import net.sf.samtools.BAMIndexMetaData;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import org.apache.commons.io.IOUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import edu.buffalo.fusim.gtf.GTFParseException;


public class BackgroundSelector implements GeneSelector {
    private static Log logger = LogFactory.getLog(BackgroundSelector.class);
    
    private List<TranscriptRecord> transcripts;
    private File geneModelFile;
    private GeneModelParser parser;

    private ArrayBlockingQueue<TranscriptRecord> queue;
    private File backgroundFile;
    private double rpkmCutoff;
    private int threads;

    public BackgroundSelector(File backgroundFile, double rpkmCutoff, int threads) {
        this.queue = new ArrayBlockingQueue<TranscriptRecord>(1000);
        this.backgroundFile = backgroundFile;
        this.rpkmCutoff = rpkmCutoff;
        this.threads = threads;
    }

    public List<TranscriptRecord> select() {
        if(transcripts == null) this.processBackgroundReads();
        return transcripts;
    }

    public List<TranscriptRecord> select(String[] filter) {
        if(filter == null || filter.length == 0) return this.select();

        if(transcripts == null) this.processBackgroundReads();

        Map<String,Boolean> filterMap = new HashMap<String,Boolean>();
        for(String f : filter) {
            filterMap.put(f,true);
        }

        List<TranscriptRecord> filteredList = new ArrayList<TranscriptRecord>();
        for(TranscriptRecord r : transcripts) {
            if(filterMap.containsKey(r.getGeneId()) 
               || filterMap.containsKey(r.getChrom())
               || filterMap.containsKey(r.getTranscriptId())) {
                filteredList.add(r);
            }
        }
        return filteredList;
    }

    protected void processBackgroundReads() {
        logger.info("Processing background reads...");
        logger.info("Computing RPKM values using " + threads + " threads...");
        long tstart = System.currentTimeMillis();
        this.transcripts = new ArrayList<TranscriptRecord>();

        GeneModelProducer producer = new GeneModelProducer();
        producer.start();

        // Leave one thread for the GTF producer
        if (threads > 1)
            threads--;
        ArrayList<GeneModelConsumer> consumers = new ArrayList<GeneModelConsumer>();

        for (int i = 0; i < threads; i++) {
            GeneModelConsumer consumer = new GeneModelConsumer(backgroundFile);
            consumer.start();
            consumers.add(consumer);
        }

        try {
            producer.join();
            for (GeneModelConsumer c : consumers) {
                c.join();
            }
        } catch (InterruptedException e) {}
        
        // Reduce
        for(GeneModelConsumer c : consumers) {
            for(TranscriptRecord t : c.getTranscripts()) {
                transcripts.add(t);
            }
            c.clearTranscripts();
        }

        long tend = System.currentTimeMillis();
        double totalTime = ((tend - tstart)/1000);
        logger.info("Finished processing background file in: "+totalTime + "s");
    }

    protected class GeneModelProducer extends Thread {
        private Log log = LogFactory.getLog(GeneModelProducer.class);

        public void run() {
            try {
                BufferedReader reader = IOUtils.toBufferedReader(new InputStreamReader(new FileInputStream(geneModelFile), "UTF-8"));

                String line = null;
                while ((line = reader.readLine()) != null) {
                    if (line.startsWith("#")) continue;
                    TranscriptRecord feature;
                    try {
                        feature = parser.parseLine(line);
                        if(feature == null) continue;

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
                log.fatal("I/O error while reading gene model file: " + e.getMessage());
            }
        }
    }

    protected class GeneModelConsumer extends Thread {
        private Log log = LogFactory.getLog(GeneModelProducer.class);
        private SAMFileReader sam;
        private int totalMappedReads;
        private List<TranscriptRecord> list = new ArrayList<TranscriptRecord>();

        public GeneModelConsumer(File bamFile) {
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
        
        public List<TranscriptRecord> getTranscripts() {
            return this.list;
        }

        public void clearTranscripts() {
            this.list = new ArrayList<TranscriptRecord>();
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
                if(rpkm > rpkmCutoff) {
                    transcript.setRPKM(rpkm);
                    list.add(transcript);
                }
            }
        }
    }

    public File getGeneModelFile() {
        return this.geneModelFile;
    }

    public void setGeneModelFile(File geneModelFile) {
        this.geneModelFile = geneModelFile;
    }

    public GeneModelParser getGeneModelParser() {
        return this.parser;
    }

    public void setGeneModelParser(GeneModelParser parser) {
        this.parser = parser;
    }

    public static void main(String[] args) throws Exception {
        BackgroundSelector b = new BackgroundSelector(new File("../sam/hg19-align.bam"), 0.2, 8);
        b.setGeneModelFile(new File("../fusim-data/refFlat.txt"));
        b.setGeneModelParser(new UCSCRefFlatParser());
        //List<TranscriptRecord> genes = b.select(new String[]{"NR_046207"});
        List<TranscriptRecord> genes = b.select();
        for(TranscriptRecord tr : genes) {
            System.out.println(tr.getTranscriptId()+": "+tr.getRPKM());
        }
    }
}
