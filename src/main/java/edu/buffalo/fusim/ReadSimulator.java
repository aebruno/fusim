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
    private Log logger = LogFactory.getLog(ReadSimulator.class);

    public static void main(String[] args) {
        ReadSimulator s = new ReadSimulator();
        s.run();
    }

    public void run() {
        // art_illumina -i fusion.txt -o testsim -l 75 -f 10 -p -m 400 -s 10 -sam
        CommandLine cmdLine = new CommandLine("art_illumina");
        cmdLine.addArgument("-i");
        cmdLine.addArgument("${file}");
        cmdLine.addArgument("-o");
        cmdLine.addArgument("testsim");
        cmdLine.addArgument("-l");
        cmdLine.addArgument("75");
        cmdLine.addArgument("-f");
        cmdLine.addArgument("10");
        cmdLine.addArgument("-p");
        cmdLine.addArgument("-m");
        cmdLine.addArgument("400");
        cmdLine.addArgument("-s");
        cmdLine.addArgument("10");
        cmdLine.addArgument("-sam");

        Map map = new HashMap();
        map.put("file", new File("fusions.fa"));
        cmdLine.setSubstitutionMap(map);

        DefaultExecutor executor = new DefaultExecutor();
        executor.setExitValue(0);
        ExecuteWatchdog watchdog = new ExecuteWatchdog(60000);
        executor.setWatchdog(watchdog);

        try {
            int exitValue = executor.execute(cmdLine);
        } catch(Exception e) {
            e.printStackTrace();
            logger.fatal("Failed to execute command: "+e.getMessage());
        }
    }
}
