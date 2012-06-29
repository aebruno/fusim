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
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/**
 * Fusim
 * 
 * @author Andrew E. Bruno
 * 
 */
public class Fusim {
    private static Log logger = LogFactory.getLog(Fusim.class);
    private Options options;

    public static void main(String[] args) {
        Fusim fusim = new Fusim();

        try {
            fusim.run(args);
        } catch (Exception e) {
            logger.fatal("Execution failed: " + e.getMessage());
            e.printStackTrace();
        }
    }
    
    public void run(String[] args) throws IOException, InterruptedException {
        buildOptions();
        
        CommandLineParser parser = new PosixParser();
        CommandLine cmd = null;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            printHelpAndExit(options, e.getMessage());
        }

        if (cmd.hasOption("h")) {
            printHelpAndExit(options);
        }
        
        if(cmd.hasOption("t")) {
            if(!"txt".equalsIgnoreCase(cmd.getOptionValue("t")) && !"fasta".equalsIgnoreCase(cmd.getOptionValue("t"))) {
                printHelpAndExit(options, "Invalid output format type. Must be txt or fasta.");
            }
        }
        
        File outputFile = null;
        if(cmd.hasOption("o")) {
            outputFile = new File(cmd.getOptionValue("o"));
        }
        
        File geneModelFile = new File(cmd.getOptionValue("g"));
        if(!geneModelFile.canRead()) {
            printHelpAndExit(options, "Please provide a valid Gene Model file");
        }
        
        if("fasta".equalsIgnoreCase(cmd.getOptionValue("t")) && !cmd.hasOption("r")) {
            printHelpAndExit(options, "missing option \"-r\". You must provide a genome reference file for FASTA output");
        }
        

        File referenceFile = null;
        
        if(cmd.hasOption("r")) {
            referenceFile = new File(cmd.getOptionValue("r"));
        }
        
        if("fasta".equalsIgnoreCase(cmd.getOptionValue("t")) && !referenceFile.canRead()) {
            printHelpAndExit(options, "Please provide a valid reference file in fasta format");
        } 
        
        if("fasta".equalsIgnoreCase(cmd.getOptionValue("t"))) {
            File referenceIndexFile = new File(referenceFile.getAbsolutePath() + ".fai");
            if(!referenceIndexFile.canRead()) {
                fatalError("Missing index file. Please index your fasta file with: samtools faidx my_genome.fa");
            }
        }

        
        int nFusions = 50;
        if(cmd.hasOption("n")) {
            try {
                nFusions = Integer.parseInt(cmd.getOptionValue("n"));
            } catch(NumberFormatException e) {
                printHelpAndExit(options, "Number of fusions (n) must be a number");
            }
        }
        

        
        File bamFile = null;
        if(cmd.hasOption("b")) {
            bamFile = new File(cmd.getOptionValue("b"));
        
            if(!bamFile.canRead()) {
                printHelpAndExit(options, "Please provide a valid BAM file");
            }
        }
        
        logger.info("------------------------------------------------------------------------");
        logger.info("Running Fusim with the following settings:");
        logger.info("------------------------------------------------------------------------");
        logger.info("Input Gene Model file: "+geneModelFile.getAbsolutePath());
        if(cmd.hasOption("b")) {
            logger.info("Background BAM file: "+bamFile.getAbsolutePath());
        }
        if(cmd.hasOption("n")) {
            logger.info("Total number of generated fusions: "+nFusions);
        }
        if(cmd.hasOption("o")) {
            logger.info("Output file: "+outputFile.getAbsolutePath());
        }
        logger.info("------------------------------------------------------------------------");
        
        FusionGenerator fg = null;
        
        if(cmd.hasOption("b")) {
            fg = new BackgroundGenerator(bamFile);
        } else {
            fg = new RandomGenerator();
        }
        
        long tstart = System.currentTimeMillis();
        List<FusionGene> fusions = fg.generate(geneModelFile, nFusions);
        long tend = System.currentTimeMillis();
        
        logger.info("Done processing.");
        double totalTime = ((tend - tstart)/1000);
        logger.info("Total processing time: " + totalTime + "s");
        
        
        OutputStream ostream = null;
        if(cmd.hasOption("o")) {
            ostream = new FileOutputStream(outputFile);
        } else {
            ostream = System.out;
        }
        
        PrintWriter out = new PrintWriter(new OutputStreamWriter(ostream, "UTF-8"));
        out.println(StringUtils.join(FusionGene.getHeader(), "\t"));
        for(FusionGene f : fusions) {
            //out.println(f);
            if(cmd.hasOption("t") && "fasta".equalsIgnoreCase(cmd.getOptionValue("t"))) {
                out.println(f.genFASTA(referenceFile, cmd.hasOption("c")));
            } else {
                out.println(f.genTXT(cmd.hasOption("c")));
            }
        }
        out.flush();
    }
    
    @SuppressWarnings("static-access")
    private void buildOptions() {
        options = new Options();
        
        options.addOption(
                OptionBuilder.withLongOpt("help")
                             .withDescription("print usage info")
                             .create("h")
            );
        options.addOption(
            OptionBuilder.withLongOpt("genemodel")
                         .withDescription("Gene Model file")
                         .hasArg()
                         .isRequired()
                         .create("g")
        );
        options.addOption(
                OptionBuilder.withLongOpt("reference")
                             .withDescription("Reference genome indexed fasta, fai")
                             .hasArg()
                             .create("r")
            );
        options.addOption(
                OptionBuilder.withLongOpt("bam")
                             .withDescription("Background BAM file")
                             .hasArg()
                             .create("b")
            );
        options.addOption(
                OptionBuilder.withLongOpt("out")
                             .withDescription("Output file")
                             .hasArg()
                             .create("o")
            );
        options.addOption(
                OptionBuilder.withLongOpt("fusions")
                             .withDescription("Total number of fusions to generate")
                             .hasArg()
                             .create("n")
            );
        options.addOption(
                OptionBuilder.withLongOpt("type")
                             .withDescription("Format of output [fasta|txt]")
                             .hasArg()
                             .create("t")
            );
        options.addOption(
                OptionBuilder.withLongOpt("cds")
                             .withDescription("Only include CDS exons")
                             .create("c")
            );
    }

    public void fatalError(String message) {
        if (message != null) logger.fatal("Fatal error: " + message);

        System.exit(1);
    }

    public void printHelpAndExit(Options options, String message) {
        if (message != null)
            logger.fatal("Usage error: " + message + "\n");
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("fusim", options);
        if (message != null) {
            System.exit(1);
        } else {
            System.exit(0);
        }
    }

    public void printHelpAndExit(Options options) {
        printHelpAndExit(options, null);
    }
}
