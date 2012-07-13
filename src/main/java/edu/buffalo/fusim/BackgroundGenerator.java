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

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import cern.colt.list.IntArrayList;
import cern.jet.random.sampling.RandomSamplingAssistant;

public class BackgroundGenerator implements FusionGenerator {
    private Log logger = LogFactory.getLog(BackgroundGenerator.class);
    private GeneSelector selector;
    private GeneSelectionMethod method;
    private List<String[]> filters;

    public List<FusionGene> generate(int nFusions, int genesPerFusion) {
        List<FusionGene> fusions = new ArrayList<FusionGene>();

        List<TranscriptRecord> transcripts = selector.select();
        if(transcripts.size() == 0) return fusions;
        
        Collections.sort(transcripts, new TranscriptCompare());
        
        // First bin genes into RPKM buckets
        GeneBins geneBins = new GeneBins(nFusions);
        geneBins.fill(transcripts);

        if(GeneSelectionMethod.BINNED.equals(method)) {
            logger.info("Generating fusions using RPKM bins...");
            for(int i = 0; i < geneBins.size(); i++) {
                IntArrayList b = geneBins.getBin(i);
                b.trimToSize();
                if(b.elements().length < genesPerFusion) {
                    logger.fatal("Not enough genes in this bin to generate fusion!");
                    continue;
                }
                List<TranscriptRecord>  genes = new ArrayList<TranscriptRecord>();
                int[] sample = RandomSamplingAssistant.sampleArray(genesPerFusion, b.elements());
                for(int s = 0; s < sample.length; s++) {
                    genes.add(transcripts.get(sample[s]));

                    // Self-fusion
                    if(sample.length == 1) {
                        genes.add(transcripts.get(sample[s]));
                    }
                }
                fusions.add(new FusionGene(genes));
            }
        } else if(GeneSelectionMethod.EMPIRICAL.equals(method)) {
            logger.info("Generating fusions based on empirical background distribution...");
            Random r = new Random();
            int[] distribution = new int[transcripts.size()];
            int index = 0;
            for(int i = 0; i < geneBins.size(); i++) {
                IntArrayList b = geneBins.getBin(i);
                b.trimToSize();
                for(int j = 0; j < b.size(); j++) {
                    distribution[index] = i;    
                    index++;
                }
            }

            for(int x = 0; x < nFusions; x++) {
                List<TranscriptRecord> genes = new ArrayList<TranscriptRecord>();

                for(int j = 0; j < genesPerFusion; j++) {
                    int binIndex = distribution[r.nextInt(distribution.length)];
                    IntArrayList b = geneBins.getBin(binIndex);
                    TranscriptRecord tr = transcripts.get(b.get(r.nextInt(b.size())));
                    genes.add(tr);

                    // Self-fusion
                    if(genesPerFusion == 1) {
                        genes.add(tr);
                    }
                }
                fusions.add(new FusionGene(genes));
            }
        } else  {
            logger.info("Generating fusions based on uniform distribution...");
            RandomGenerator rg = new RandomGenerator();
            rg.setGeneSelector(selector);
            rg.setFilters(filters);
            rg.setGeneSelectionMethod(method);
            fusions = rg.generate(nFusions, genesPerFusion);
        }
        
        return fusions;
    }
  
    protected class TranscriptCompare implements Comparator<TranscriptRecord> {
        public int compare(TranscriptRecord o1, TranscriptRecord o2) {
            return Double.compare(o1.getRPKM(), o2.getRPKM());
        }
    }
    
    protected class GeneBins {
        private IntArrayList[] bins;
        private int binsize;
        
        public GeneBins(int binsize) {
            this.binsize = binsize;
        }
        
        public void fill(List<TranscriptRecord> values) {
            if(values.size() < binsize) {
                this.binsize = 1;
            }

            this.bins = new IntArrayList[binsize];
            for(int i = 0; i < binsize; i++) {
                bins[i] = new IntArrayList();
            }

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
