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

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * GENCODE GTF Parser
 * 
 * See: http://www.gencodegenes.org/gencodeformat.html
 * 
 * GTF File Format: http://mblab.wustl.edu/GTF22.html
 * 
 * @author Andrew E. Bruno
 * 
 */
public class GTFParser {
    private static final Pattern ATTRIBUTE_PATTERN = Pattern
            .compile("^\\s*(.+)\\s(.+)$");

    public Feature parseLine(String line) throws GTFParseException {
        if (line == null || line.startsWith("#"))
            return null;

        GencodeFeature record = new GencodeFeature();

        String[] fields = line.split("\t");
        record.seqname = fields[0];
        record.source = fields[1];
        record.featureType = FeatureType.fromString(fields[2]);

        try {
            record.start = Integer.valueOf(fields[3]);
        } catch (NumberFormatException e) {
            throw new GTFParseException("Invalid integer value for start", e);
        }

        try {
            record.end = Integer.valueOf(fields[4]);
        } catch (NumberFormatException e) {
            throw new GTFParseException("Invalid integer value for end", e);
        }

        record.strand = Strand.fromString(fields[6]);

        try {
            record.score = Double.valueOf(fields[5]);
            record.frame = Integer.valueOf(fields[7]);
        } catch (NumberFormatException ignored) {
        }

        record.attributes = new HashMap<String, String>();
        
        if (fields.length >= 8 && fields[8] != null) {
            String[] attrs = fields[8].split(";");

            for (String variableString : attrs) {
                Matcher m = ATTRIBUTE_PATTERN.matcher(variableString);
                if (m.matches()) {
                    String key = m.group(1).trim();
                    String val = m.group(2).trim();
                    val = val.replaceAll("\"", "");
                    if (val.length() > 0) {
                        record.attributes.put(key, val);
                    }
                }
            }
        }
        
        //XXX should probably throw an exception if these are null
        record.geneId = record.getAttribute("gene_id");
        record.transcriptId = record.getAttribute("transcript_id");
        
        //XXX These fields are specific to Gencode. Consider doing more type checking
        record.geneType = record.getAttribute("gene_type");
        record.geneStatus = record.getAttribute("gene_status");
        record.geneName = record.getAttribute("gene_name");
        record.transcriptType = record.getAttribute("transcript_type");
        record.transcriptStatus = record.getAttribute("transcript_status");
        record.transcriptName = record.getAttribute("transcript_name");
        try {
            record.level = Integer.valueOf(record.getAttribute("level"));
        } catch(NumberFormatException ignored) { }
        

        return record;
    }

    //XXX test only. remove soon
    public static void main(String[] args) throws Exception {
        GTFParser parser = new GTFParser();

        FileInputStream in = new FileInputStream("test.gtf");
        BufferedReader reader = new BufferedReader(new InputStreamReader(in));
        String line = null;
        while ((line = reader.readLine()) != null) {
            Feature f = parser.parseLine(line);
            System.out.print(f);
        }
    }

}
