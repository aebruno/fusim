package edu.buffalo.fusim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import edu.buffalo.fusim.gtf.Feature;
import edu.buffalo.fusim.gtf.FeatureType;
import edu.buffalo.fusim.gtf.GTFParseException;
import edu.buffalo.fusim.gtf.GTFParser;
import edu.buffalo.fusim.gtf.Strand;

public class GTF2RefFlat {
    private static Log logger = LogFactory.getLog(GTF2RefFlat.class);
    
    private GTFParser parser = new GTFParser();
    private Map<String, List<String>> map;
    
    
    /**
     * This is a port of the mkFromGroupedGxf(..) function from genePred.c from 
     * Kent source utilities
     * http://genomewiki.cse.ucsc.edu/index.php/Kent_source_utilities
     */
    public void convert(File gtfFile, File outFile) throws IOException {
        logger.info("Converting GTF File: "+gtfFile.getAbsolutePath());
        long tstart = System.currentTimeMillis();
        
        map = new HashMap<String, List<String>>();
        buildMap(gtfFile);
        PrintWriter output = new PrintWriter(new OutputStreamWriter(new FileOutputStream(outFile), "UTF-8"));
        
        for(String id : map.keySet()) {
            TranscriptData data = null;
            for(String line : map.get(id)) {
                Feature feature;
                try {
                    feature = parser.parseLine(line);
                } catch (GTFParseException e) {
                    logger.fatal("Invalid GTF format. Could not parse line: "+e.getMessage());
                    continue;
                }
                if(data == null) {
                    data = new TranscriptData(feature.getTranscriptId(), feature.getGeneId(), feature.getSeqname(), feature.getStrand());
                }
                
                data.addFeature(feature);
            }

            Collections.sort(data.getFeatures(), new FeatureCompare());

            int exonCount = 0;
            int stopCodonStart = -1; 
            int stopCodonEnd = -1;
            int cdsStart = Integer.MAX_VALUE;
            int cdsEnd = Integer.MIN_VALUE;
            int txStart = Integer.MAX_VALUE;
            int txEnd = Integer.MIN_VALUE;
            boolean haveStartCodon = false;
            boolean haveStopCodon = false;

            for(Feature f : data.getFeatures()) {
                if (txStart > f.getStart()) txStart = f.getStart();
                if (txEnd < f.getEnd()) txEnd = f.getEnd();

                if(FeatureType.EXON.equals(f.getFeatureType())) {
                    exonCount++;
                }
                if(FeatureType.CDS.equals(f.getFeatureType())) {
                    if (f.getStart() < cdsStart)
                        cdsStart = f.getStart();
                    if (f.getEnd() > cdsEnd)
                        cdsEnd = f.getEnd();
                }
                if(FeatureType.START_CODON.equals(f.getFeatureType())) 
                    haveStartCodon = true;
                if(FeatureType.STOP_CODON.equals(f.getFeatureType()))  {
                    /* stop_codon can be split, need bounds for adjusting CDS below */
                    if ((stopCodonStart < 0) || (f.getStart() < stopCodonStart))
                        stopCodonStart = f.getStart();
                    if ((stopCodonEnd < 0) || (f.getEnd() > stopCodonEnd))
                        stopCodonEnd = f.getEnd();

                    haveStopCodon = true;
                }
            }

            if (exonCount == 0) continue;

            if (cdsStart > cdsEnd) {
                /* no cds annotated */
                cdsStart = 0;
                cdsEnd = 0;
            } else if (stopCodonStart >= 0) {
                /* adjust CDS to include stop codon as in GTF */
                if (Strand.FORWARD.equals(data.getStrand())) {
                    if (stopCodonEnd > cdsEnd) cdsEnd = stopCodonEnd;
                } else {
                    if (stopCodonStart < cdsStart) cdsStart = stopCodonStart;
                }
            }

            if(cdsStart > cdsEnd) {
                cdsStart = txStart;
                cdsEnd = txEnd;
            }

            /* adjust tx range to include stop codon */
            if (Strand.FORWARD.equals(data.getStrand()) && (txEnd == stopCodonStart))
                 txEnd = stopCodonEnd;
             else if (Strand.REVERSE.equals(data.getStrand()) && (txStart == stopCodonEnd))
                 txStart = stopCodonStart;
            
            int[] exonStarts = new int[exonCount];
            int[] exonEnds = new int[exonCount];

            int i = -1; /* before first exon */
            /* fill in exons, merging overlaping and adjacent exons */
            for(Feature f : data.getFeatures()) {
                if(FeatureType.EXON.equals(f.getFeatureType()) || FeatureType.CDS.equals(f.getFeatureType())) {
                    if ((i < 0) || (f.getStart() > exonEnds[i])) {
                        /* start a new exon */
                        ++i;
                        assert(i < exonCount);
                        exonStarts[i] = f.getStart();
                        exonEnds[i] = f.getEnd();
                    } else {
                        /* overlap, extend exon, picking the largest of ends */
                        assert(i < exonCount);
                        assert(f.getStart() >= exonStarts[i]);
                        if (f.getEnd() > exonEnds[i])
                            exonEnds[i] = f.getEnd();
                    }
                }
            }

            exonCount = i+1;

            StringBuffer buf = new StringBuffer();
            
            buf.append(
            StringUtils.join(new String[]{
                    data.getGeneId(), 
                    data.getTranscriptId(), 
                    data.getChrom(), 
                    data.getStrand().toString(),
                    ""+txStart,
                    ""+txEnd,
                    ""+cdsStart,
                    ""+cdsEnd,
                    ""+exonCount,
                    StringUtils.join(ArrayUtils.toObject(exonStarts), ","),
                    StringUtils.join(ArrayUtils.toObject(exonEnds), ",")
            }, "\t"));
            
            output.println(buf.toString());
        }

        output.flush();
        
        long tend = System.currentTimeMillis();
        double totalTime = ((tend - tstart)/1000);
        logger.info("Finished conversion: "+totalTime + "s");
        logger.info("Output written to: "+outFile.getAbsolutePath());
    }
    
    protected class FeatureCompare implements Comparator<Feature> {
        public int compare(Feature o1, Feature o2) {
            return Double.compare(o1.getStart(), o2.getStart());
        }
    }
    
    private void buildMap(File gtfFile) throws IOException {
        try {
            BufferedReader reader = IOUtils.toBufferedReader(new InputStreamReader(new FileInputStream(gtfFile), "UTF-8"));
            
            String line = null;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#")) continue;
                Feature feature = parser.parseLine(line);
                
                List<String> lines = map.get(feature.getTranscriptId());
                if(lines == null) {
                    lines = new ArrayList<String>();
                }
                
                lines.add(line);
                map.put(feature.getTranscriptId(), lines);
            }
        } catch (IOException e) {
            throw new RuntimeException("Failed to read gene modle file", e);
        } catch (GTFParseException e) {
            throw new RuntimeException("Failed to parse gene model file", e);
        }
    }

    //XXX test only. remove soon
    public static void main(String[] args) throws Exception {
        GTF2RefFlat g = new GTF2RefFlat();
        
        g.convert(new File("../fusim-data/gencodev11.gtf"), new File("tmp.rg"));
    }
    
    
    private class TranscriptData {
        private List<Feature> features = new ArrayList<Feature>();
        private String transcriptId;
        private String geneId;
        private String chrom;
        private Strand strand;
        
        public TranscriptData(String transcriptId, String geneId, String chrom, Strand strand) {
            this.transcriptId = transcriptId;
            this.geneId = geneId;
            this.chrom = chrom;
            this.strand = strand;
        }
        
        public void addFeature(Feature f) {
            features.add(f);
        }

        public List<Feature> getFeatures() {
            return features;
        }

        public String getTranscriptId() {
            return transcriptId;
        }

        public String getGeneId() {
            return geneId;
        }

        public String getChrom() {
            return chrom;
        }

        public Strand getStrand() {
            return strand;
        }
        
        
    }
}
