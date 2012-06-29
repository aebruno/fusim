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
import java.io.FileInputStream;
import java.io.InputStreamReader;

import cern.colt.Arrays;

public class GeneModelParser {
    public TranscriptRecord parseLine(String line) {
        if (line == null || line.startsWith("#"))
            return null;

        String[] fields = line.split("\t");
        TranscriptRecord record = TranscriptRecord.fromGeneModel(fields);

        return record;
    }

    //XXX test only. remove soon
    public static void main(String[] args) throws Exception {
        GeneModelParser parser = new GeneModelParser();

        FileInputStream in = new FileInputStream("data/refGene.txt");
        BufferedReader reader = new BufferedReader(new InputStreamReader(in));
        String line = null;
        while ((line = reader.readLine()) != null) {
            TranscriptRecord f = parser.parseLine(line);
            //if(f!=null && f.getCodingExons().size() == 0) System.out.println(f.getTranscriptId());
            if(f != null && ("NR_002206".equals(f.getTranscriptId()) ||   "NR_026911".equals(f.getTranscriptId()))) {
                System.out.print(f);
                System.out.println(Arrays.toString(f.generateExonBreak(false, true)));
            } 
        }
    }
}
