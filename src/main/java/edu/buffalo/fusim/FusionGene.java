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

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.collections.Transformer;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;

import edu.buffalo.fusim.gtf.Strand;

public class FusionGene {

    private List<TranscriptRecord> genes = new ArrayList<TranscriptRecord>();
    private List<FusionOption> options = new ArrayList<FusionOption>();
    private FusionType fusionType = FusionType.HYBRID;
    private String geneId;
    private String transcriptId;

    public FusionGene(List<TranscriptRecord> transcripts) { 
        for(TranscriptRecord tr : transcripts) {
            genes.add(tr);
        }
        this.setIds();
    }
    
    public FusionGene(TranscriptRecord gene1, TranscriptRecord gene2) {
        this.genes.add(gene1);
        this.genes.add(gene2);
        this.setIds();
    }

    public FusionGene(TranscriptRecord gene1, TranscriptRecord gene2, TranscriptRecord gene3) {
        this.genes.add(gene1);
        this.genes.add(gene2);
        this.genes.add(gene3);
        this.setIds();
    }
    
    public String outputFasta(List<int []> breaks, File reference, boolean cdsExonsOnly, boolean fixOrientation) {
        return this.outputFasta(breaks, reference, cdsExonsOnly, fixOrientation, 0);
    }

    public String outputFasta(List<int []> breaks, File reference, boolean cdsExonsOnly, boolean fixOrientation, int foreignInsertionLen) {
        ExtractSeq extractSeq = new ExtractSeq(reference);

        StringBuffer fasta = new StringBuffer();
        fasta.append(">ref|"+this.getTranscriptId()
                     +" fusionGene="+this.getGeneId()
                     +" fusionType="+this.getFusionType()
                     +" fusionOptions="+StringUtils.join(this.options, ","));

        List<StringBuffer> seqs = new ArrayList<StringBuffer>();
        for(int b = 0; b < breaks.size(); b++) {
            int[] exons = breaks.get(b); 
            TranscriptRecord gene = genes.get(b);
            int breakno = b+1;
            fasta.append(" chrom"+breakno+"="+gene.getChrom()
                         +" strand"+breakno+"="+gene.getStrand()
                         +" exonIndex"+breakno+"="+StringUtils.join(ArrayUtils.toObject(exons), ",")
                         //+" exons"+breakno+"="
                         );
        
            StringBuffer breakSeq = new StringBuffer();
            for(int i = 0; i < exons.length; i++) {
                int[] exon = gene.getExons(cdsExonsOnly).get(exons[i]);
                //fasta.append((exon[0]+1)+"-"+exon[1]);
                //if(i != (exons.length-1)) fasta.append(",");
                breakSeq.append(extractSeq.fetch(gene.getChrom(), gene.getStrand(), exon[0]+1, exon[1]));
            }
            seqs.add(breakSeq);
        }

        fasta.append("\n");
        
        Strand normStrand = genes.get(0).getStrand();
        
        for(int i = 0; i < genes.size(); i++) {
            if(fixOrientation && !normStrand.equals(genes.get(i).getStrand())) {
                fasta.append(ExtractSeq.reverseComplement(seqs.get(i)));
            } else {
                fasta.append(seqs.get(i).toString());   
            }

            if(foreignInsertionLen > 0) {
                if(i == 0 || (i == 1 && genes.size() == 3)) {
                    fasta.append(ExtractSeq.randomSequence(foreignInsertionLen));
                }
            }
        }

        
        return fasta.toString();
    }
    
    private List<String> createTxtColumns(TranscriptRecord gene, int[] breaks, boolean cdsExonsOnly) {
        List<String> cols = new ArrayList<String>();
        cols.add(this.getGeneId());
        cols.add(gene.getGeneId());
        cols.add(gene.getTranscriptId());
        cols.add(gene.getChrom());
        cols.add(gene.getStrand().toString());
        cols.add(""+breaks.length);

        int[] exonStarts = new int[breaks.length];
        int[] exonEnds = new int[breaks.length];
        int exonBases = 0;
        for(int i = 0; i < breaks.length; i++) {
            int[] exon = gene.getExons(cdsExonsOnly).get(breaks[i]);
            exonStarts[i] = (exon[0]+1);
            exonEnds[i] = exon[1];
            exonBases += exon[1]-exon[0];
        }

        cols.add(""+exonBases);
        cols.add(StringUtils.join(ArrayUtils.toObject(breaks), ","));
        cols.add(StringUtils.join(ArrayUtils.toObject(exonStarts), ","));
        cols.add(StringUtils.join(ArrayUtils.toObject(exonEnds), ","));
        cols.add(this.fusionType.toString());
        cols.add(StringUtils.join(this.options, ","));

        return cols;
    }

    public String outputText(List<int[]> breaks, boolean cdsExonsOnly) {
        StringBuffer txt = new StringBuffer();

        for(int i = 0; i < breaks.size(); i++) {
            int[] exons = breaks.get(i);
            TranscriptRecord gene = genes.get(i);
            txt.append(StringUtils.join(this.createTxtColumns(gene, exons, cdsExonsOnly), "\t")+"\n");
        }

        return txt.toString();
    }

    private void setIds() {
        this.geneId = StringUtils.join(
                            CollectionUtils.collect(
                                this.genes, 
                                new Transformer() {
                                    public Object transform(Object g) {
                                        return ((TranscriptRecord)g).getGeneId();
                                    }
                                }
                            ), "-");

        this.transcriptId = StringUtils.join(
                            CollectionUtils.collect(
                                this.genes, 
                                new Transformer() {
                                    public Object transform(Object g) {
                                        return ((TranscriptRecord)g).getTranscriptId();
                                    }
                                }
                            ), "-");
    }

    public String getGeneId() {
        return this.geneId;
    }

    public String getTranscriptId() {
        return this.transcriptId;
    }
    
    public static String[] getHeader() {
        return new String[]{
                "fusionGene", "geneName", "name", "chrom", "strand", "exonCount",
                "exonBases", "exonIndexes", "exonStarts", "exonEnds", "fusionType", "fusionOptions"
                };
    }
    
    public String toString() {
        String str = "";
        for(int i = 0; i < genes.size(); i++) {
            TranscriptRecord gene = genes.get(i);
            int geneno = i+1;
            str += "---- FUSION "+geneno+" -------\n";
            str += gene.toString();
        }
        return str;
    }

    public TranscriptRecord getGene(int index) {
        return genes.get(index);
    }

    public int size() {
        return genes.size();
    }

    public void addGene(TranscriptRecord gene) {
        genes.add(gene);
    }

    public FusionType getFusionType() {
        return this.fusionType;
    }

    public void setFusionType(FusionType fusionType) {
        this.fusionType = fusionType;
    }

    public List<FusionOption> getFusionOptions() {
        return this.options;
    }

    public void addOption(FusionOption option) {
        this.options.add(option);
    }
}
