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
        
        CommandLineParser clParser = new PosixParser();
        CommandLine cmd = null;
        try {
            cmd = clParser.parse(options, args);
        } catch (ParseException e) {
            printHelpAndExit(options, e.getMessage());
        }

        if (cmd.hasOption("h")) {
            printHelpAndExit(options);
        }

        if(cmd.hasOption("z")) {
            if(!cmd.hasOption("i")) {
                printHelpAndExit(options, "Please specify a path to a GTF/GFF file for conversion with option -i");
            }
            if(!cmd.hasOption("o")) {
                printHelpAndExit(options, "Please specify an output filename with option -o");
            }

            File outFile = new File(cmd.getOptionValue("o"));
            File gtfFile = new File(cmd.getOptionValue("i"));
            if(!gtfFile.canRead()) {
                printHelpAndExit(options, "Can't read input GTF file");
            }

            GTF2RefFlat gtf2Flat = new GTF2RefFlat();
            
            gtf2Flat.convert(gtfFile, outFile);
            System.exit(0);
        }
        
        if(!cmd.hasOption("g")) {
            printHelpAndExit(options, "Please specify a path to a gene model file with option -g");
        }
        
        File geneModelFile = new File(cmd.getOptionValue("g"));
        if(!geneModelFile.canRead()) {
            printHelpAndExit(options, "Can't read Gene Model file");
        }
        
        PrintWriter textOutput = null;
        if(cmd.hasOption("t")) {
            if("-".equals(cmd.getOptionValue("t"))) {
                textOutput = new PrintWriter(new OutputStreamWriter(System.out, "UTF-8"));
            } else {
                textOutput = new PrintWriter(new OutputStreamWriter(new FileOutputStream(cmd.getOptionValue("t")), "UTF-8"));
            }
        } 
        
        PrintWriter fastaOutput = null;
        if(cmd.hasOption("f")) {
            if("-".equals(cmd.getOptionValue("f"))) {
                fastaOutput = new PrintWriter(new OutputStreamWriter(System.out, "UTF-8"));
            } else {
                fastaOutput = new PrintWriter(new OutputStreamWriter(new FileOutputStream(cmd.getOptionValue("f")), "UTF-8"));
            }
        }
        
        // Default to TXT output
        if(fastaOutput == null && textOutput == null) {
            textOutput = new PrintWriter(new OutputStreamWriter(System.out, "UTF-8"));
        }

        if(cmd.hasOption("u") && !cmd.hasOption("f")) {
            printHelpAndExit(options, "You must provide an output FASTA file for simulating Illumina reads");
        }
        
        if(cmd.hasOption("f") && !cmd.hasOption("r")) {
            printHelpAndExit(options, "You must provide an indexed (.fai) genome reference file for FASTA output using option \"-r\".");
        }

        File referenceFile = null;
        
        if(cmd.hasOption("r")) {
            referenceFile = new File(cmd.getOptionValue("r"));
        }
        
        if(cmd.hasOption("f") && !referenceFile.canRead()) {
            printHelpAndExit(options, "Please provide a valid reference file in fasta format");
        } 
        
        if(cmd.hasOption("f")) {
            File referenceIndexFile = new File(referenceFile.getAbsolutePath() + ".fai");
            if(!referenceIndexFile.canRead()) {
                fatalError("Missing index file. Please index your fasta file with: samtools faidx my_genome.fa");
            }
        }

        
        int nFusions = 10;
        if(cmd.hasOption("n")) {
            try {
                nFusions = Integer.parseInt(cmd.getOptionValue("n"));
            } catch(NumberFormatException e) {
                printHelpAndExit(options, "Number of fusions (-n) must be a number");
            }
        }
        
        int nReadThrough = 0;
        if(cmd.hasOption("x")) {
            try {
                nReadThrough = Integer.parseInt(cmd.getOptionValue("x"));
            } catch(NumberFormatException e) {
                printHelpAndExit(options, "Number of read through fusion genes (-x) must be a number");
            }
        }

        int nThreads = Runtime.getRuntime().availableProcessors();
        if(cmd.hasOption("p")) {
            try {
                nThreads = Integer.parseInt(cmd.getOptionValue("p"));
            } catch(NumberFormatException e) {
                printHelpAndExit(options, "Number of threads to spawn (-p) must be a number");
            }
        }

        int readLength = 75;
        if(cmd.hasOption("l")) {
            try {
                readLength = Integer.parseInt(cmd.getOptionValue("l"));
            } catch(NumberFormatException e) {
                printHelpAndExit(options, "Read length (-l) must be a number");
            }
        }

        int meanFrag = 400;
        if(cmd.hasOption("m")) {
            try {
                meanFrag = Integer.parseInt(cmd.getOptionValue("m"));
            } catch(NumberFormatException e) {
                printHelpAndExit(options, "Mean DNA fragment length (-m) must be a number");
            }
        }

        int readCoverage = 10;
        if(cmd.hasOption("d")) {
            try {
                readCoverage = Integer.parseInt(cmd.getOptionValue("d"));
            } catch(NumberFormatException e) {
                printHelpAndExit(options, "Read coverage (-l) must be a number");
            }
        }

        String artPrefix = "fusion-reads";
        if(cmd.hasOption("y")) {
            artPrefix = cmd.getOptionValue("y");
        }

        String artPath = ReadSimulator.DEFAULT_ART_BIN;
        if(cmd.hasOption("a")) {
            artPath = cmd.getOptionValue("a");
        }

        double rpkmCutoff = 0.2;
        if(cmd.hasOption("k")) {
            try {
                rpkmCutoff = Double.parseDouble(cmd.getOptionValue("k"));
                if(rpkmCutoff < 0 || rpkmCutoff > 1) throw new NumberFormatException();
            } catch(NumberFormatException e) {
                printHelpAndExit(options, "RPKM cutoff (-k) must be 0 < cutoff < 1");
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
        if(cmd.hasOption("n")) {
            logger.info("Total number of generated fusions: "+nFusions);
        }
        if(cmd.hasOption("x")) {
            logger.info("Total number of read through genes: "+nReadThrough);
        }
        if(cmd.hasOption("t")) {
            logger.info("Text Output file: "+cmd.getOptionValue("t"));
        }
        if(cmd.hasOption("f")) {
            logger.info("Fasta Output file: "+cmd.getOptionValue("f"));
        }
        if(cmd.hasOption("b")) {
            logger.info("-- Generating fusions based on background dataset --");
            logger.info("Background BAM file: "+bamFile.getAbsolutePath());
            logger.info("RPKM cutoff: "+rpkmCutoff);
            logger.info("Number of threads: "+nThreads);
        }
        if(cmd.hasOption("u")) {
            logger.info("-- Simulating Illumina reads using ART --");
            logger.info("ART Path: "+artPath);
            logger.info("ART output prefix: "+artPrefix);
            logger.info("Read length: "+readLength);
            logger.info("Mean DNA Fragment length: "+meanFrag);
            logger.info("The fold of read coverage: "+readCoverage);
            logger.info("Paired-end reads?: "+(cmd.hasOption("e") ? "Yes" : "No"));
        }
        logger.info("------------------------------------------------------------------------");
        
        GeneModelParser parser = new UCSCRefFlatParser();
        FusionGenerator fg = null;
        
        if(cmd.hasOption("b")) {
            fg = new BackgroundGenerator(bamFile, parser, rpkmCutoff, nThreads);
        } else {
            fg = new RandomGenerator(parser);
        }
        
        logger.info("Starting fusion gene simulation...");
        long tstart = System.currentTimeMillis();
        List<FusionGene> fusions = fg.generate(geneModelFile, nFusions);
        long tend = System.currentTimeMillis();
        
        logger.info("Simulation complete.");
        double totalTime = ((tend - tstart)/1000);
        logger.info("Total processing time: " + totalTime + "s");
        
        // Generate any read through fusion genes
        if(nReadThrough > 0) {
            logger.info("Generating read through genes...");
            ReadThroughGenerator rt = new ReadThroughGenerator(parser);
            List<FusionGene> rtFusions = rt.generate(geneModelFile, nReadThrough);
            fusions.addAll(rtFusions);
        }
        
        if(textOutput != null) {
            textOutput.println(StringUtils.join(FusionGene.getHeader(), "\t"));
        }
        
        for(FusionGene f : fusions) {
            //out.println(f);
            
            // First half of gene 1
            int[] break1 = f.getGene1().generateExonBreak(true, cmd.hasOption("c"));
            
            // Second half of gene2
            int[] break2 = f.getGene2().generateExonBreak(false, cmd.hasOption("c"));
            
            if(textOutput != null) {
                textOutput.println(f.genTXT(break1, break2, cmd.hasOption("c")));
            }
            
            if(fastaOutput != null) {
                fastaOutput.println(f.genFASTA(break1, break2, referenceFile, cmd.hasOption("c")));
            }
        }
        
 
        
        if(textOutput != null) textOutput.flush();
        if(fastaOutput != null) fastaOutput.flush();

        if(cmd.hasOption("u")) {
            logger.info("Simulating Illumina reads using ART...");
            ReadSimulator s = new ReadSimulator();
            s.run(artPath, new File(cmd.getOptionValue("f")), artPrefix, readLength, meanFrag, readCoverage, cmd.hasOption("e"));
        }

        logger.info("Fusim run complete. Goodbye!");
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
                OptionBuilder.withLongOpt("fusions")
                             .withDescription("Total number of fusions to generate")
                             .hasArg()
                             .create("n")
            );
        options.addOption(
                OptionBuilder.withLongOpt("text")
                             .withDescription("File name of text output")
                             .hasArg()
                             .create("t")
            );
        options.addOption(
                OptionBuilder.withLongOpt("fasta")
                             .withDescription("File name of fasta output")
                             .hasArg()
                             .create("f")
            );
        options.addOption(
                OptionBuilder.withLongOpt("cds")
                             .withDescription("Only include CDS exons")
                             .create("c")
            );
        options.addOption(
                OptionBuilder.withLongOpt("readthrough")
                             .withDescription("Number of read through fusion genes")
                             .hasArg()
                             .create("x")
            );
        options.addOption(
                OptionBuilder.withLongOpt("cutoff")
                             .withDescription("RPKM cutoff when using background BAM file")
                             .hasArg()
                             .create("k")
            );
        options.addOption(
                OptionBuilder.withLongOpt("threads")
                             .withDescription("Number of threads to spawn when processing background BAM file")
                             .hasArg()
                             .create("p")
            );
        options.addOption(
                OptionBuilder.withLongOpt("convert")
                             .withDescription("Convert GTF/GFF to refFlat (genePred) format")
                             .create("z")
            );
        options.addOption(
                OptionBuilder.withLongOpt("illumina")
                             .withDescription("Simulate Illumina Reads with ART")
                             .create("u")
            );
        options.addOption(
                OptionBuilder.withLongOpt("gtf")
                             .withDescription("Input GTF file for conversion")
                             .hasArg()
                             .create("i")
            );
        options.addOption(
                OptionBuilder.withLongOpt("output")
                             .withDescription("Output refFlat file for conversion")
                             .hasArg()
                             .create("o")
            );
        options.addOption(
                OptionBuilder.withLongOpt("art")
                             .withDescription("Path to ART binary for simulating Illumina reads from fusion genes")
                             .hasArg()
                             .create("a")
            );
        options.addOption(
                OptionBuilder.withLongOpt("prefix")
                             .withDescription("Prefix of output files from ART for simulating Illumina reads")
                             .hasArg()
                             .create("y")
            );
        options.addOption(
                OptionBuilder.withLongOpt("readlength")
                             .withDescription("Length of reads to be simulated from ART for simulating Illumina reads")
                             .hasArg()
                             .create("l")
            );
        options.addOption(
                OptionBuilder.withLongOpt("meanfrag")
                             .withDescription("Mean size of DNA fragments for paired-end reads from ART for simulating Illumina reads")
                             .hasArg()
                             .create("m")
            );
        options.addOption(
                OptionBuilder.withLongOpt("coverage")
                             .withDescription("The fold of read coverage to be simulated using ART for simulating Illumina reads")
                             .hasArg()
                             .create("d")
            );
        options.addOption(
                OptionBuilder.withLongOpt("paired")
                             .withDescription("Simulate paired-end reads using ART for simulating Illumina reads")
                             .create("e")
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
