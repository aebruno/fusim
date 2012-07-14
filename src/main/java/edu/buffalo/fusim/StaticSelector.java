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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.HashMap;
import java.util.Map;

import edu.buffalo.fusim.gtf.GTFParseException;

import org.apache.commons.io.IOUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

public class StaticSelector implements GeneSelector {
    private static Log logger = LogFactory.getLog(StaticSelector.class);
    
    private List<TranscriptRecord> transcripts;
    private File geneModelFile;
    private GeneModelParser parser;

    public List<TranscriptRecord> select() {
        if(transcripts == null) this.parseTranscripts();
        return transcripts;
    }

    public List<TranscriptRecord> select(String[] filter) {
        if(filter == null || filter.length == 0) return this.select();

        if(transcripts == null) this.parseTranscripts();

        Map<String,Boolean> filterMap = new HashMap<String,Boolean>();
        for(String f : filter) {
            filterMap.put(f,true);
        }

        List<TranscriptRecord> filteredList = new ArrayList<TranscriptRecord>();
        for(TranscriptRecord r : transcripts) {
            if(filterMap.containsKey(r.getGeneId()) 
                || filterMap.containsKey(r.getChrom())
                || filterMap.containsKey(r.getTranscriptId())) {
                filteredList.add(r);
            }
        }

        if(filteredList.size() == 0) {
            throw new RuntimeException("No transcripts found using filter: "+Arrays.toString(filter));
        }
        return filteredList;
    }

    private void parseTranscripts() {
        logger.info("Parsing gene model file...");
        long tstart = System.currentTimeMillis();
        this.transcripts = new ArrayList<TranscriptRecord>();
        
        try {
            BufferedReader reader = IOUtils.toBufferedReader(new InputStreamReader(new FileInputStream(geneModelFile), "UTF-8"));
            
            String line = null;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#")) continue;
                
                TranscriptRecord record = parser.parseLine(line);
                if(record == null) continue;
                
                transcripts.add(record);
            }
        } catch (IOException e) {
            throw new RuntimeException("Failed to read gene modle file", e);
        } catch (GTFParseException e) {
            throw new RuntimeException("Failed to parse gene model file", e);
        }

        long tend = System.currentTimeMillis();
        double totalTime = ((tend - tstart)/1000);
        logger.info("Finished parsing gene model file in: "+totalTime + "s");
        if(transcripts.size() == 0) {
            throw new RuntimeException("No transcripts found! Can't generate fusions without transcripts!");
        }
    }

    public File getGeneModelFile() {
        return this.geneModelFile;
    }

    public void setGeneModelFile(File geneModelFile) {
        this.geneModelFile = geneModelFile;
    }

    public GeneModelParser getGeneModelParser() {
        return this.parser;
    }

    public void setGeneModelParser(GeneModelParser parser) {
        this.parser = parser;
    }

    public static void main(String[] args) throws Exception {
        StaticSelector s = new StaticSelector();
        s.setGeneModelFile(new File("../fusim-data/refFlat.txt"));
        s.setGeneModelParser(new UCSCRefFlatParser());
        List<TranscriptRecord> genes = s.select(new String[]{"NR_046207"});
        for(TranscriptRecord tr : genes) {
            System.out.println(tr.getTranscriptId());
        }
    }
}
