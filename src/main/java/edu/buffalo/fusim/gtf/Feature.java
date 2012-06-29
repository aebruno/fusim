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

import java.util.Map;


/**
 * Feature class encapsulating a "feature" or "annotation" from a GTF/GFF file
 * The format is defined here:
 * 
 * http://mblab.wustl.edu/GTF22.html
 * 
 * Descriptions of fields were copied directly from above URL
 * 
 * @author Andrew E. Bruno
 * 
 */
public class Feature {
    protected String seqname;
    protected String source;
    protected FeatureType featureType;
    protected int start;
    protected int end;
    protected Double score;
    protected Strand strand;
    protected Integer frame;
    protected Map<String, String> attributes;
    protected String geneId;
    protected String transcriptId;
    
    public Feature() {

    }

    /**
     * The name of the sequence. Commonly, this is the chromosome ID or contig
     * ID. Note that the coordinates used must be unique within each sequence
     * name in all GTFs for an annotation set.
     */
    public String getSeqname() {
        return seqname;
    }

    /**
     * Genomic strand {+,-}
     */
    public Strand getStrand() {
        return strand;
    }

    /**
     * Genomic feature type.
     */
    public FeatureType getFeatureType() {
        return featureType;
    }

    /**
     * start coordinate of the feature relative to the beginning of the
     * sequence. Numbering is 1 based.
     * 
     * @return
     */
    public int getStart() {
        return start;
    }

    /**
     * end coordinate of the feature. Numbering is 1 based.
     */
    public int getEnd() {
        return end;
    }

    /**
     * a unique label indicating where the annotations came from. typically the
     * name of either a prediction program or a public database
     */
    public String getSource() {
        return source;
    }

    /**
     * a degree of confidence in the feature's existence and coordinates. The
     * value of this field has no global scale but may have relative
     * significance when the source field indicates the prediction program used
     * to create this annotation. It may be a floating point number or integer.
     */
    public Double getScore() {
        return score;
    }

    /**
     * 0 indicates that the feature begins with a whole codon at the 5' most
     * base. 1 means that there is one extra base (the third base of a codon)
     * before the first whole codon and 2 means that there are two extra bases
     * (the second and third bases of the codon) before the first codon. Note
     * that for reverse strand features, the 5' most base is the end coordinate.
     */
    public Integer getFrame() {
        return frame;
    }

    /**
     * Additional attributes for genomic feature.
     */
    public Map<String, String> getAttributes() {
        return attributes;
    }
    
    /**
     * Returns the value of a specific attribute
     * @param key attribute key
     * @return attribute value
     */
    public String getAttribute(String key) {
        return attributes.get(key);
    }

    /**
     * A globally unique identifier for the genomic locus of the transcript. If
     * empty, no gene is associated with this feature
     */
    public String getGeneId() {
        return geneId;
    }

    /**
     * A globally unique identifier for the predicted transcript. If empty, no
     * transcript is associated with this feature
     */
    public String getTranscriptId() {
        return transcriptId;
    }

    public String toString() {
        String str = "[\n"; 
        str += "Seqname: "+seqname + ':' + start + '-' + end + " " + strand + "\n";
        str += "Source: " +source + "\n";
        str += "Feature Type: " + featureType + "\n";
        str += "Score: "+ score + "\n";
        str += "Frame: "+frame + "\n";
        
        if(attributes != null) {
            str += "Attributes: \n";
            for(String key : attributes.keySet()) {
                str += "    "+key+" = "+attributes.get(key)+"\n";
            }
        }
        
        str += "]\n";
        
        return str;
    }
}
