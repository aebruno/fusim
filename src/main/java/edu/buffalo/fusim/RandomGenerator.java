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
import java.util.List;
import java.util.Random;

public class RandomGenerator implements FusionGenerator {
    
    private GeneSelector selector;
    private GeneSelectionMethod method;
    private List<String[]> filters;

    public List<FusionGene> generate(int nFusions, int genesPerFusion) {
        List<FusionGene> list = new ArrayList<FusionGene>();
        Random r = new Random();

        for(int n = 0; n < nFusions; n++) {
            List<TranscriptRecord> genes = new ArrayList<TranscriptRecord>();
            for(int i = 0; i < genesPerFusion; i++) {
                List<TranscriptRecord> transcripts = null;
                if(filters != null && i <= filters.size()-1) {
                    transcripts = selector.select(filters.get(i));
                } else {
                    transcripts = selector.select();
                }
                TranscriptRecord t = transcripts.get(r.nextInt(transcripts.size()));
                genes.add(t);

                // self-fusion
                if(genesPerFusion == 1) {
                    genes.add(t);
                }
            }
            list.add(new FusionGene(genes));
        }

        return list;
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
