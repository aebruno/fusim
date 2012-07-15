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
import java.util.Map;
import java.util.HashMap;
import java.util.Random;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import cern.colt.list.IntArrayList;
import cern.jet.random.sampling.RandomSamplingAssistant;

public class IntraChromGenerator implements FusionGenerator {
    private Log logger = LogFactory.getLog(IntraChromGenerator.class);
    private GeneSelector selector;
    private GeneSelectionMethod method;
    private List<String[]> filters;

    private static String[] chroms = new String[]{
        "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
        "chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
        "chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"
    };

    public List<FusionGene> generate(int nFusions, int genesPerFusion) {
        List<FusionGene> fusions = new ArrayList<FusionGene>();

        Random r = new Random();

        RandomGenerator rg = new RandomGenerator();
        rg.setGeneSelector(selector);
        rg.setGeneSelectionMethod(method);

        for(int i = 0; i < nFusions; i++) {
            String chr = chroms[r.nextInt(chroms.length)];
            List<String[]> list = new ArrayList<String[]>();
            list.add(new String[]{chr});
            list.add(new String[]{chr});
            rg.setFilters(list);
            fusions.addAll(rg.generate(1, genesPerFusion));
        }
        
        return fusions;
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
