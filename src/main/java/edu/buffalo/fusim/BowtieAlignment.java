/*
 * Copyright 2012 Andrew E. Bruno <aebruno2@buffalo.edu> *
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
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.exec.CommandLine;
import org.apache.commons.exec.DefaultExecutor;
import org.apache.commons.exec.ExecuteWatchdog;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

public class BowtieAlignment {
    public static final String DEFAULT_BOWTIE_BIN = "bowtie";

    private Log logger = LogFactory.getLog(BowtieAlignment.class);

    public static void main(String[] args) {
        BowtieAlignment b = new BowtieAlignment();
        //b.run(BowtieAlignment.DEFAULT_BOWTIE_BIN, "hg19", new File("fusion-reads1.fq"), new File("fusion-reads2.fq"), new File("fusion.sam"), 8);
        b.run(BowtieAlignment.DEFAULT_BOWTIE_BIN, "hg19", new File("fusion-reads1.fq"), new File("fusion.sam"), 8);
    }

    // Single-end 
    public void run(String bowtiePath, String index, File reads, File outFile, int threads) {
        this.run(bowtiePath, index, reads, null, outFile, threads);
    }

    // Paried-end
    public void run(String bowtiePath, String index, File read1, File read2, File outFile, int threads) {
        // Example bowtie call:
        //   bowtie -t -v 2 -p 12 -m 10 -S
        CommandLine cmdLine = new CommandLine(bowtiePath);

        // print wall-clock time taken by search phases
        cmdLine.addArgument("-t");
        // eport end-to-end hits w/ <=v mismatches; ignore qualities
        cmdLine.addArgument("-v");
        cmdLine.addArgument("2");
        // number of alignment threads to launch
        cmdLine.addArgument("-p");
        cmdLine.addArgument(""+threads);
        // suppress all alignments if > <int> exist
        cmdLine.addArgument("-m");
        cmdLine.addArgument("10");
        // write hits in SAM format
        cmdLine.addArgument("-S");
        // bowtie index
        cmdLine.addArgument("${index}");
        if(read2 != null) {
            // fastq1
            cmdLine.addArgument("-1");
            cmdLine.addArgument("${read1}");
            // fastq2
            cmdLine.addArgument("-2");
            cmdLine.addArgument("${read2}");
        } else {
            cmdLine.addArgument("${read1}");
        }
        // output SAM file
        cmdLine.addArgument("${outFile}");

        Map map = new HashMap();
        map.put("index", index);
        map.put("read1", read1);
        if(read2 != null) {
            map.put("read2", read2);
        }
        map.put("outFile", outFile);
        cmdLine.setSubstitutionMap(map);

        DefaultExecutor executor = new DefaultExecutor();
        executor.setExitValue(0);
        // Never timeout 
        ExecuteWatchdog watchdog = new ExecuteWatchdog(ExecuteWatchdog.INFINITE_TIMEOUT);
        executor.setWatchdog(watchdog);

        try {
            int exitValue = executor.execute(cmdLine);
        } catch(Exception e) {
            logger.fatal("Failed to execute bowtie for aligning simulated Illumina reads from ART: "+e.getMessage());
        }
    }
}
