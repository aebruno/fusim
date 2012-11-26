package edu.buffalo.fusim;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

import edu.buffalo.fusim.gtf.GTFParseException;

public class ReadThroughGenerator implements FusionGenerator {

    private GeneSelector selector;
    private GeneSelectionMethod method;
    private List<String[]> filters;

    public List<FusionGene> generate(int nFusions, int genesPerFusion) {
        List<FusionGene> list = new ArrayList<FusionGene>();

        //XXX ignoring filters for now..
        List<TranscriptRecord> transcripts = selector.select();
        if(transcripts.size() < genesPerFusion) return list;

        Collections.sort(transcripts, new TranscriptCompare());
        
        Random r = new Random();

        //XXX we only support 2 genes per ReadThrough fusion
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
            
            return Integer.valueOf(o1.getTxStart()).compareTo(Integer.valueOf(o2.getTxStart()));
        }
    }

    public void setGeneSelector(GeneSelector selector) {
        this.selector = selector;
    }

    public GeneSelector getGeneSelector() {
        return this.selector;
    }

    public void setGeneSelectionMethod(GeneSelectionMethod method) {
        this.method = method;
    }

    public GeneSelectionMethod getGeneSelectionMethod() {
        return this.method;
    }

    public void setFilters(List<String[]> filters) {
        this.filters = filters;
    }

    public List<String[]> getFilters() {
        return this.filters;
    }
}
