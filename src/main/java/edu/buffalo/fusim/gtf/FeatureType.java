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
 * Enum for encapsulating the various Genomic feature types found in the Gencode
 * data. See http://www.gencodegenes.org/data.html
 * 
 * This is rather strict. For example, if a GTF file doesn't contain one of the 
 * types listed in this class it will throw a parse error. To support any type
 * of GTF file consider just removing this enum and using strings.
 * 
 * @author Andrew E. Bruno
 * 
 */
public enum FeatureType {
    CDS, EXON, GENE, START_CODON, STOP_CODON, TRANSCRIPT, UTR, UTR5, UTR3, SELENOCYSTEINE, INTER, INTER_CNS, INTRON_CNS;

    public static FeatureType fromString(String str) throws GTFParseException {
        if (str.equalsIgnoreCase("cds")) {
            return FeatureType.CDS;
        } else if (str.equalsIgnoreCase("exon")) {
            return FeatureType.EXON;
        } else if (str.equalsIgnoreCase("gene")) {
            return FeatureType.GENE;
        } else if (str.equalsIgnoreCase("Selenocysteine")) {
            return FeatureType.SELENOCYSTEINE;
        } else if (str.equalsIgnoreCase("start_codon")) {
            return FeatureType.START_CODON;
        } else if (str.equalsIgnoreCase("stop_codon")) {
            return FeatureType.STOP_CODON;
        } else if (str.equalsIgnoreCase("transcript")) {
            return FeatureType.TRANSCRIPT;
        } else if (str.equalsIgnoreCase("UTR")) {
            return FeatureType.UTR;
        } else if (str.equalsIgnoreCase("5utr")) {
            return FeatureType.UTR5;
        } else if (str.equalsIgnoreCase("3utr")) {
            return FeatureType.UTR3;
        } else if (str.equalsIgnoreCase("inter")) {
            return FeatureType.INTER;
        } else if (str.equalsIgnoreCase("inter_cns")) {
            return FeatureType.INTER_CNS;
        } else if (str.equalsIgnoreCase("intron_cns")) {
            return FeatureType.INTRON_CNS;
        } else {
            throw new GTFParseException("Invalid feature type '" + str + "'.");
        }
    }
}
