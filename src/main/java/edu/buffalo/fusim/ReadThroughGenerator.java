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

public class ReadThroughGenerator implements FusionGenerator {
    private GeneModelParser parser;

    public ReadThroughGenerator() {
        this.parser = new GeneModelParser();
    }
        
    public List<FusionGene> generate(File gtfFile, int nFusions) {
        List<TranscriptRecord> transcripts = this.parseTranscripts(gtfFile);
        Collections.sort(transcripts, new TranscriptCompare());
        
        Random r = new Random();

        List<FusionGene> list = new ArrayList<FusionGene>();
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
        
        // Sorty by Chrom then txStart
        public int compare(TranscriptRecord o1, TranscriptRecord o2) {
            int res = o1.getChrom().compareTo(o2.getChrom());
            if (res != 0) return res;
            
            return Double.compare(o1.getTxStart(), o2.getTxStart());
        }
    }
    
    private List<TranscriptRecord> parseTranscripts(File gtfFile) {
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
                transcripts.add(record);
            }
        } catch (IOException e) {
            throw new RuntimeException("Failed to read gene modle file", e);
        }
        
        return transcripts;
    }
    
    public static void main(String[] args) throws Exception {
        ReadThroughGenerator g = new ReadThroughGenerator();
        System.out.println(g.generate(new File("../fusim-data/refGene.txt"), 1));
    }

}
