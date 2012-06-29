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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;

import edu.buffalo.fusim.gtf.Feature;
import edu.buffalo.fusim.gtf.FeatureType;
import edu.buffalo.fusim.gtf.GTFParseException;
import edu.buffalo.fusim.gtf.GTFParser;
import edu.buffalo.fusim.gtf.Strand;

public class GTF2RefFlat {
    
    private GTFParser parser = new GTFParser();
    private Map<String, TranscriptData> map;
    
    
    public void convert(File gtfFile, File outFile) throws IOException {
        map = new HashMap<String, TranscriptData>();
        buildMap(gtfFile);
        
        for(String id : map.keySet()) {
            TranscriptData data = map.get(id);
            
            // Skip transcripts without tx, cds, or exon data
            if(data.getTx().size() == 0 && data.getCds().size() == 0) continue;
            if(data.getExon().size() == 0) continue;
            
            Feature tx = null;
            Feature cds = null;
            if(data.getTx().size() > 0) {
                tx = data.getTx().get(0);
            }
            if(data.getCds().size() > 0) {
                cds = data.getCds().get(0);
            }
            
            if(tx == null) tx = cds;
            if(cds == null) cds = tx;
            
            Collections.sort(data.getExon(), new ExonCompare());

            int[] exonStarts = new int[data.getExon().size()];
            int[] exonEnds = new int[data.getExon().size()];
            for(int i = 0; i < data.getExon().size(); i++) {
                exonStarts[i] = data.getExon().get(i).getStart();
                exonEnds[i] = data.getExon().get(i).getEnd();
            }
            
            StringBuffer buf = new StringBuffer();
            
            buf.append(
            StringUtils.join(new String[]{
                    data.getGeneId(), 
                    data.getTranscriptId(), 
                    data.getChrom(), 
                    data.getStrand().toString(),
                    ""+tx.getStart(),
                    ""+tx.getEnd(),
                    ""+cds.getStart(),
                    ""+cds.getEnd(),
                    ""+data.getExon().size(),
                    StringUtils.join(ArrayUtils.toObject(exonStarts), ","),
                    StringUtils.join(ArrayUtils.toObject(exonEnds), ",")
            }, "\t"));
            
            System.out.println(buf.toString());
        }
    }
    
    protected class ExonCompare implements Comparator<Feature> {
        public int compare(Feature o1, Feature o2) {
            return Double.compare(o1.getStart(), o2.getStart());
        }
    }
    
    private void buildMap(File gtfFile) throws IOException {
        // XXX Use Java NIO for reading File. Requires Java 7
        try {
            BufferedReader reader = Files.newBufferedReader(FileSystems
                    .getDefault().getPath(gtfFile.getAbsolutePath()),
                    Charset.forName("UTF-8"));
            
            String line = null;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#")) continue;
                Feature feature = parser.parseLine(line);
                
                TranscriptData data = map.get(feature.getTranscriptId());
                if(data == null) {
                    data = new TranscriptData(feature.getTranscriptId(), feature.getGeneId(), feature.getSeqname(), feature.getStrand());
                }
                
                if(FeatureType.EXON.equals(feature.getFeatureType())) {
                    data.addExon(feature);
                } else if(FeatureType.CDS.equals(feature.getFeatureType())) {
                    data.addCds(feature);
                } else if(FeatureType.TRANSCRIPT.equals(feature.getFeatureType())) {
                    data.addTx(feature);
                }
                
                map.put(feature.getTranscriptId(), data);
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
        private List<Feature> cds = new ArrayList<Feature>();
        private List<Feature> exon = new ArrayList<Feature>();
        private List<Feature> tx = new ArrayList<Feature>();
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
        
        public void addCds(Feature f) {
            cds.add(f);
        }

        public void addExon(Feature f) {
            exon.add(f);
        }

        public void addTx(Feature f) {
            tx.add(f);
        }

        public List<Feature> getCds() {
            return cds;
        }

        public List<Feature> getExon() {
            return exon;
        }

        public List<Feature> getTx() {
            return tx;
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
