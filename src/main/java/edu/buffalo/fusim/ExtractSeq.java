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
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.util.SequenceUtil;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import edu.buffalo.fusim.gtf.Strand;

public class ExtractSeq {
    private static Log logger = LogFactory.getLog(ExtractSeq.class);

    private static String[] bases = new String[]{"A", "C", "T", "G"};
    
    private static Map<Character,Character> symbolMap = new HashMap<Character,Character>();
    static {
        symbolMap.put('a', 't');
        symbolMap.put('A', 'T');
        symbolMap.put('c', 'g');
        symbolMap.put('C', 'G');
        symbolMap.put('g', 'c');
        symbolMap.put('G', 'C');
        symbolMap.put('t', 'a');
        symbolMap.put('T', 'A');
    }
    
    private ReferenceSequenceFile ref;
    
    public ExtractSeq(File path) {
        this.ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(path); 
    }
    
    public String fetch(String chrom, Strand strand, int start, int stop) {
        //logger.info("Fetching sequence contig: "+chrom+":"+start+"-"+stop+" "+strand);
        StringBuilder buff = new StringBuilder();
        
        ReferenceSequence seq = ref.getSubsequenceAt(chrom, start, stop);
        byte[] bases = seq.getBases();
        if (Strand.REVERSE.equals(strand)) SequenceUtil.reverseComplement(bases);
            
        for (int i=0; i<bases.length; ++i) {
            buff.append((char)bases[i]);
        }  
        
        return buff.toString();
    }
    
    public static StringBuffer reverseComplement(StringBuffer seq) {
        StringBuffer revc = new StringBuffer();
        
        seq.reverse();
        for(int i = 0; i < seq.length(); i++) {
            revc.append(symbolMap.get(seq.charAt(i))); 
        }
        
        return revc;
    }

    public static StringBuffer randomSequence(int maxLen) {
        Random r = new Random();

        // Chose a random length
        int length = r.nextInt(maxLen)+1;

        StringBuffer seq = new StringBuffer();
        for(int i = 0; i < length; i++) {
            seq.append(bases[r.nextInt(bases.length)]);
        }

        return seq;
    }
    
    public static void main(String[] args) {
        //ExtractSeq s = new ExtractSeq(new File("data/hg19.fa"));
        //System.out.println(s.fetch("chr17", Strand.FORWARD, 76210398, 76210508));
        
        //StringBuffer buf = new StringBuffer("ACTG");
        //System.out.println(ExtractSeq.reverseComplement(buf).toString());
        
        System.out.println(ExtractSeq.randomSequence(10));
    }

}
