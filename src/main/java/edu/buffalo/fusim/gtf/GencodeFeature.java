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

package edu.buffalo.fusim.gtf;


/**
 * Feature class specific to GENCODE
 * 
 * The format is defined here:
 * 
 * http://www.gencodegenes.org/gencodeformat.html
 * 
 * 
 * @author Andrew E. Bruno
 * 
 */
public class GencodeFeature extends Feature {

    protected String geneType;
    protected String geneStatus;
    protected String geneName;
    protected String transcriptType;
    protected String transcriptStatus;
    protected String transcriptName;
    protected Integer level;

    public GencodeFeature() {
    }
    
    /**
     * Chromosome
     */
    public String getChrom() {
        return seqname;
    }

    /**
     * Gene Biotype
     */
    public String getGeneType() {
        return geneType;
    }

    /**
     * Gene status: {KNOWN, NOVEL, PUTATIVE}
     */
    public String getGeneStatus() {
        return geneStatus;
    }

    /**
     * Gene Name
     */
    public String getGeneName() {
        return geneName;
    }

    /**
     * Transcript Biotype
     */
    public String getTranscriptType() {
        return transcriptType;
    }

    /**
     * Transcript status {KNOWN, NOVEL, PUTATIVE}
     */
    public String getTranscriptStatus() {
        return transcriptStatus;
    }

    /**
     * Transcript name
     */
    public String getTranscriptName() {
        return transcriptName;
    }

    /**
     * 1 (verified loci), 2 (manually annotated loci), 3 (automatically
     * annotated loci)
     */
    public Integer getLevel() {
        return level;
    }

}
