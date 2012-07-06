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

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;

public class FusionGene {
    
    private TranscriptRecord gene1;
    private TranscriptRecord gene2;
    
    public FusionGene(TranscriptRecord gene1, TranscriptRecord gene2) {
        this.gene1 = gene1;
        this.gene2 = gene2;
    }
    
    public String outputFasta(int[] break1, int[] break2, File reference, boolean cdsExonsOnly, boolean fixOrientation) {
        ExtractSeq extractSeq = new ExtractSeq(reference);

        StringBuffer fasta = new StringBuffer();
        fasta.append(">ref|"+gene1.getTranscriptId()+"-"+gene2.getTranscriptId()
                     +" fusionGene="+gene1.getGeneId()+"-"+gene2.getGeneId()
                     +" exons1="+StringUtils.join(ArrayUtils.toObject(break1), ",")
                     +" break1="+gene1.getChrom()+":");
        
        StringBuffer break1seq = new StringBuffer();
        for(int i = 0; i < break1.length; i++) {
            int[] exon = gene1.getExons(cdsExonsOnly).get(break1[i]);
            fasta.append((exon[0]+1)+"-"+exon[1]);
            if(i != (break1.length-1)) fasta.append(",");
            break1seq.append(extractSeq.fetch(gene1.getChrom(), gene1.getStrand(), exon[0]+1, exon[1]));
        }
        
        fasta.append(" strand1="+gene1.getStrand());
        fasta.append(" exons2="+StringUtils.join(ArrayUtils.toObject(break2), ","));
        fasta.append(" break2="+gene2.getChrom()+":");
        
        StringBuffer break2seq = new StringBuffer();
        for(int i = 0; i < break2.length; i++) {
            int[] exon = gene2.getExons(cdsExonsOnly).get(break2[i]);
            fasta.append((exon[0]+1)+"-"+exon[1]);
            if(i != (break2.length-1)) fasta.append(",");
            break2seq.append(extractSeq.fetch(gene2.getChrom(), gene2.getStrand(), exon[0]+1, exon[1]));
        }
        
        fasta.append(" strand2="+gene2.getStrand()+"\n");
        fasta.append(break1seq.toString());
        
        if(!gene1.getStrand().equals(gene2.getStrand())) {
            if(fixOrientation) {
                fasta.append(ExtractSeq.reverseComplement(break2seq));
            } else {
                fasta.append(break2seq.toString());
            }
        } else {
            fasta.append(break2seq.toString());   
        }

        
        return fasta.toString();
    }
    
    private List<String> createTxtColumns(TranscriptRecord gene, int[] breaks, boolean cdsExonsOnly) {
        List<String> cols = new ArrayList<String>();
        cols.add(gene1.getGeneId()+"-"+gene2.getGeneId());
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

        return cols;
    }

    public String outputText(int[] break1, int[] break2, boolean cdsExonsOnly) {
        StringBuffer txt = new StringBuffer();

        txt.append(StringUtils.join(this.createTxtColumns(gene1, break1, cdsExonsOnly), "\t")+"\n");
        txt.append(StringUtils.join(this.createTxtColumns(gene2, break2, cdsExonsOnly), "\t")+"\n");

        return txt.toString();
    }
    
    public static String[] getHeader() {
        return new String[]{
                "fusionGene", "geneName", "name", "chrom", "strand", "exonCount",
                "exonBases", "exonIndexes", "exonStarts", "exonEnds"
                };
    }
    
    public String toString() {
        String str = "---- FUSION 1 -------\n";
        str += gene1.toString();
        str += "---- FUSION 2 -------\n";
        str += gene2.toString();
        return str;
    }

    public TranscriptRecord getGene1() {
        return gene1;
    }

    public TranscriptRecord getGene2() {
        return gene2;
    }
   
}
