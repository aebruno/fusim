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

public class ReadSimulator {
    public static final String DEFAULT_ART_BIN = "art_illumina";

    private Log logger = LogFactory.getLog(ReadSimulator.class);

    public static void main(String[] args) {
        ReadSimulator s = new ReadSimulator();
        s.run(ReadSimulator.DEFAULT_ART_BIN, new File("fusions.fa"), "test-fusion-sims", 75, 400, 10, true);
    }

    public void run(String artBinPath, File fusionFile, String outputPrefix, int readLength, int meanFragSize, int readCoverage, boolean pairedEnd) {
        // Example art call:
        //   art_illumina -i fusion.txt -o testsim -l 75 -f 10 -p -m 400 -s 10 
        CommandLine cmdLine = new CommandLine(artBinPath);

        // the filename of input DNA reference
        cmdLine.addArgument("-i");
        cmdLine.addArgument("${file}");
        // the prefix of output files
        cmdLine.addArgument("-o");
        cmdLine.addArgument("${outputPrefix}");
        // the length of reads to be simulated
        cmdLine.addArgument("-l");
        cmdLine.addArgument(""+readLength);
        // the fold of read coverage to be simulated
        cmdLine.addArgument("-f");
        cmdLine.addArgument(""+readCoverage);
        if(pairedEnd) {
            // indicate a paired-end read simulation
            cmdLine.addArgument("-p");
            // the mean size of DNA fragments for paired-end simulations
            cmdLine.addArgument("-m");
            cmdLine.addArgument(""+meanFragSize);
            // the standard deviation of DNA fragment size for paired-end simulations.
            cmdLine.addArgument("-s");
            cmdLine.addArgument("10");
        }
        // quite - turn off end of run summary
        cmdLine.addArgument("-q");

        Map map = new HashMap();
        map.put("file", fusionFile);
        map.put("outputPrefix", outputPrefix);
        cmdLine.setSubstitutionMap(map);

        DefaultExecutor executor = new DefaultExecutor();
        executor.setExitValue(0);
        // Timeout after 5 minutes
        ExecuteWatchdog watchdog = new ExecuteWatchdog(300000);
        executor.setWatchdog(watchdog);

        try {
            int exitValue = executor.execute(cmdLine);
        } catch(Exception e) {
            logger.fatal("Failed to execute ART for simulating Illumina reads: "+e.getMessage());
        }
    }
}
