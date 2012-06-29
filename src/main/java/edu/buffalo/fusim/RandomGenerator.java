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
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import edu.buffalo.fusim.gtf.GTFParseException;

public class RandomGenerator implements FusionGenerator {
    private GeneModelParser parser;

    public RandomGenerator(GeneModelParser parser) {
        this.parser = parser;
    }

    public List<FusionGene> generate(File gtfFile, int nFusions) {
        List<FusionGene> list = new ArrayList<FusionGene>();
        Random r = new Random();
        
        try {
            RandomAccessFile raf = new RandomAccessFile(gtfFile, "r");
            for (int i = 0; i < nFusions; i++) {
                list.add(new FusionGene(this.getRandomTranscript(r, raf, (int)gtfFile.length()),
                                        this.getRandomTranscript(r, raf, (int)gtfFile.length())));
            }
        } catch (IOException e) {
            throw new RuntimeException("Failed to generate fusion genes: "+e.getMessage());
        }

        return list;
    }

    private TranscriptRecord getRandomTranscript(Random r,
            RandomAccessFile raf, int length) throws IOException {
        int randomByte = r.nextInt(length);
        raf.seek((long) randomByte);
        raf.readLine();

        TranscriptRecord record = null;
        for (int tries = 0; tries < 5; tries++) {
            String line = raf.readLine();
            try {
                record = parser.parseLine(line);
            } catch (GTFParseException e) {
                // XXX ignored for now;
            }
            if (record != null)
                break;
        }

        if (record == null) {
            throw new IOException(
                    "Failed to find random transcript after 5 attempts");
        }

        return record;
    }

    public static void main(String[] args) throws Exception {
        RandomGenerator g = new RandomGenerator(new UCSCRefFlatParser());
        System.out.println(g.generate(new File("data/refGene.txt"), 1));
    }
}
