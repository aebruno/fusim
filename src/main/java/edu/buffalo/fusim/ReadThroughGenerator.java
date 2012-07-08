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
import java.util.Map;
import java.util.Random;

import edu.buffalo.fusim.gtf.GTFParseException;

public class ReadThroughGenerator implements FusionGenerator {
    private GeneModelParser parser;

    public ReadThroughGenerator(GeneModelParser parser) {
        this.parser = parser;
    }
        
    public List<FusionGene> generate(File gtfFile, int nFusions, GeneSelectionMethod method, Map<String, Boolean> limit) {
        List<FusionGene> list = new ArrayList<FusionGene>();

        List<TranscriptRecord> transcripts = this.parseTranscripts(gtfFile, limit);
        if(transcripts.size() < 2) {
            return list;
        }
        Collections.sort(transcripts, new TranscriptCompare());
        
        Random r = new Random();

        for (int i = 0; i < nFusions; i++) {
            int index = r.nextInt(transcripts.size()-1);
            
            TranscriptRecord gene1 = transcripts.get(index);
            TranscriptRecord gene2 = null;
            
            for(int j = index+1; j < transcripts.size(); j++) {
                TranscriptRecord t = transcripts.get(j);
                if(!gene1.getGeneId().equals(t.getGeneId())) {
                    gene2 = t;
                    break;
                }
            }
            
            if(gene2 == null) {
                for(int j = index-1; j >= 0; j--) {
                    TranscriptRecord t = transcripts.get(j);
                    if(!gene1.getGeneId().equals(t.getGeneId())) {
                        gene2 = gene1;
                        gene2 = t;
                        break;
                    }
                }
            }
            
            list.add(new FusionGene(gene1, gene2));
        }
        
        return list;
    }
    
    protected class TranscriptCompare implements Comparator<TranscriptRecord> {
        
        // Sort by Chrom then txStart
        public int compare(TranscriptRecord o1, TranscriptRecord o2) {
            int res = o1.getChrom().compareTo(o2.getChrom());
            if (res != 0) return res;
            
            return Integer.compare(o1.getTxStart(), o2.getTxStart());
        }
    }

    protected static class TranscriptLenCompare implements Comparator<TranscriptRecord> {
        
        public int compare(TranscriptRecord o1, TranscriptRecord o2) {
            return Integer.compare(o1.getExonBases(), o2.getExonBases());
        }
    }
    
    protected List<TranscriptRecord> parseTranscripts(File gtfFile, Map<String, Boolean> limit) {
        List<TranscriptRecord> transcripts = new ArrayList<TranscriptRecord>();
        
        // XXX Use Java NIO for reading File. Requires Java 7
        try {
            BufferedReader reader = Files.newBufferedReader(FileSystems
                    .getDefault().getPath(gtfFile.getAbsolutePath()),
                    Charset.forName("UTF-8"));
            
            String line = null;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#")) continue;
                
                TranscriptRecord record = parser.parseLine(line);
                if(record == null) continue;
                
                //XXX skip the haplotypes and unassembled chroms
                if(record.getChrom().contains("_")) continue;

                if(limit != null && 
                  !limit.containsKey(record.getGeneId()) &&
                  !limit.containsKey(record.getTranscriptId())) continue;

                transcripts.add(record);
            }
        } catch (IOException e) {
            throw new RuntimeException("Failed to read gene modle file", e);
        } catch (GTFParseException e) {
            throw new RuntimeException("Failed to parse gene model file", e);
        }
        
        return transcripts;
    }
    
    public static void main(String[] args) throws Exception {
        ReadThroughGenerator g = new ReadThroughGenerator(new UCSCRefFlatParser());
        //System.out.println(g.generate(new File("../fusim-data/refGene.txt"), 1, GeneSelectionMethod.UNIFORM, null));
        List<TranscriptRecord> transcripts = g.parseTranscripts(new File("../fusim-data/refFlat.txt"), null);
        Collections.sort(transcripts, new TranscriptLenCompare());
        for(int i = 0; i < 100; i++) {
            System.out.println(transcripts.get(i));
        }
    }

}
