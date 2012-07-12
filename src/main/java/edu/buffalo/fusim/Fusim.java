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
import java.io.InputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;

import net.sf.picard.sam.MergeSamFiles;

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

        if (cmd.hasOption("v")) {
            printVersionAndExit();
        }

        if (cmd.hasOption("h") || cmd.getOptions().length == 0) {
            printHelpAndExit(options);
        }

        if(cmd.hasOption("z")) {
            if(!cmd.hasOption("i")) {
                printHelpAndExit(options, "Please specify a path to a GTF/GFF file for conversion with option -i");
            }
            if(!cmd.hasOption("o")) {
                printHelpAndExit(options, "Please specify an output filename with option -O");
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

        double rpkmCutoff = 0.2;
        if(cmd.hasOption("k")) {
            try {
                rpkmCutoff = Double.parseDouble(cmd.getOptionValue("k"));
                if(rpkmCutoff < 0 || rpkmCutoff > 1) throw new NumberFormatException();
            } catch(NumberFormatException e) {
                printHelpAndExit(options, "RPKM cutoff (-k) must be 0 < cutoff < 1");
            }
        }
        
        GeneSelectionMethod geneSelectioMethod = GeneSelectionMethod.UNIFORM;
        if(cmd.hasOption("m")) {
            GeneSelectionMethod sm = GeneSelectionMethod.fromString(cmd.getOptionValue("m"));
            if(sm == null) {
                printHelpAndExit(options, "Invalid gene selection method: "+cmd.getOptionValue("m"));
            }
            geneSelectioMethod = sm;
        }

        Map<String, Boolean> limit = null;
        if(cmd.hasOption("l")) {
            limit = new HashMap<String,Boolean>();
            String[] limits = cmd.getOptionValue("l").split(",");
            if(limits.length < 2) {
                printHelpAndExit(options, "Must provide a limit of at least two genes (ex. gene1,gene2,..): "+cmd.getOptionValue("l"));
            }
            for(int i = 0; i < limits.length; i++) {
                limit.put(limits[i], true);
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
        logger.info("Reference Gene Model: "+geneModelFile.getAbsolutePath());
        logger.info("Number of simulated fusion genes: "+nFusions);
        logger.info("Number of read through genes: "+nReadThrough);
        logger.info("Auto-correct orientation: "+(cmd.hasOption("a") ? "yes" : "no"));
        logger.info("Allow fusion genes outside of ORF: "+(cmd.hasOption("d") ? "yes" : "no"));
        logger.info("Force fusion breaks on exon boundries: "+(cmd.hasOption("e") ? "yes" : "no"));
        if(cmd.hasOption("t")) {
            logger.info("Text Output: "+("-".equals(cmd.getOptionValue("t")) ? "<stdout>" : cmd.getOptionValue("t")));
        }
        if(cmd.hasOption("f")) {
            logger.info("Fasta Output: "+("-".equals(cmd.getOptionValue("f")) ? "<stdout>" : cmd.getOptionValue("f")));
        }
        if(!cmd.hasOption("f") && !cmd.hasOption("t")) {
            logger.info("Text Output: <stdout>");
        }
        if(cmd.hasOption("b")) {
            logger.info("-- Generating fusions based on background dataset --");
            logger.info("Background BAM file: "+bamFile.getAbsolutePath());
            logger.info("RPKM cutoff: "+rpkmCutoff);
            logger.info("Number of threads: "+nThreads);
            logger.info("Gene selection method: "+geneSelectioMethod.toString());
        }
        logger.info("------------------------------------------------------------------------");
        
        GeneModelParser parser = new UCSCRefFlatParser(cmd.hasOption("e"), cmd.hasOption("c"), limit);
        FusionGenerator fg = null;
        
        if(cmd.hasOption("b")) {
            fg = new BackgroundGenerator(bamFile, parser, rpkmCutoff, nThreads, limit);
        } else {
            fg = new RandomGenerator(parser, limit);
        }
        
        logger.info("Starting fusion gene simulation...");
        long tstart = System.currentTimeMillis();
        List<FusionGene> fusions = fg.generate(geneModelFile, nFusions, geneSelectioMethod);
        long tend = System.currentTimeMillis();
        
        logger.info("Simulation complete.");
        double totalTime = ((tend - tstart)/1000);
        logger.info("Total processing time: " + totalTime + "s");
        
        // Generate any read through fusion genes
        if(nReadThrough > 0) {
            logger.info("Generating read through genes...");
            ReadThroughGenerator rt = new ReadThroughGenerator(parser);
            List<FusionGene> rtFusions = rt.generate(geneModelFile, nReadThrough, geneSelectioMethod);
            fusions.addAll(rtFusions);
        }

        if(fusions.size() == 0) {
            fatalError("No genes found to simulate fusions!");    
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

            // Keep ORF (don't allow out of frame) and allow splitting of exons
            if(!cmd.hasOption("d") && !cmd.hasOption("e")) {
                // Ensure fusion gene is within ORF
                int break1BaseCount = 0;
                for(int i = 0; i < break1.length; i++) {
                    int[] exon = f.getGene1().getExons(cmd.hasOption("c")).get(break1[i]);
                    break1BaseCount += exon[1]-exon[0];
                }

                int break2BaseCount = 0;
                for(int i = 0; i < break2.length; i++) {
                    int[] exon = f.getGene2().getExons(cmd.hasOption("c")).get(break2[i]);
                    break2BaseCount += exon[1]-exon[0];
                }

                if((break1BaseCount+break2BaseCount) % 3 != 0) {
                    // Auto adjust 
                    int break1BaseAdjust = 0;
                    int break2BaseAdjust = 0;
                    while(break1BaseCount % 3 != 0) {
                        break1BaseAdjust++;
                        break1BaseCount--;
                    }
                    while(break2BaseCount % 3 != 0) {
                        break2BaseAdjust++;
                        break2BaseCount--;
                    }

                    f.getGene1().getExons(cmd.hasOption("c")).get(break1[0])[1] -= break1BaseAdjust;
                    f.getGene2().getExons(cmd.hasOption("c")).get(break2[0])[1] -= break2BaseAdjust;

                }
            } else if(cmd.hasOption("e") && !cmd.hasOption("d")) {
                // Keep ORF (don't allow out of frame) and don't allow splitting of exons (keep exon boundries)
                // Break gene 1 on exons boundries
                break1 = f.getGene1().generateExonBoundryBreak(cmd.hasOption("c"));
                
                // Break gene 2 on exons boundries
                break2 = f.getGene2().generateExonBoundryBreak(cmd.hasOption("c"));

            }
            
            if(textOutput != null) {
                textOutput.print(f.outputText(break1, break2, cmd.hasOption("c")));
            }
            
            if(fastaOutput != null) {
                fastaOutput.println(f.outputFasta(break1, break2, referenceFile, cmd.hasOption("c"), cmd.hasOption("a")));
            }
        }
        
 
        
        if(textOutput != null) textOutput.flush();
        if(fastaOutput != null) fastaOutput.flush();

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
            OptionBuilder.withLongOpt("gene-model")
                         .withDescription("Gene Model file in refFlat format")
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
                OptionBuilder.withLongOpt("background-reads")
                             .withDescription("Path to BAM file containing background reads. Genes will be selected for fusions according to the read profile of the background reads.")
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
                OptionBuilder.withLongOpt("text-output")
                             .withDescription("File name of text output")
                             .hasArg()
                             .create("t")
            );
        options.addOption(
                OptionBuilder.withLongOpt("fasta-output")
                             .withDescription("File name of FASTA output")
                             .hasArg()
                             .create("f")
            );
        options.addOption(
                OptionBuilder.withLongOpt("cds-only")
                             .withDescription("Only include CDS exons")
                             .create("c")
            );
        options.addOption(
                OptionBuilder.withLongOpt("read-through")
                             .withDescription("Number of read through fusion genes")
                             .hasArg()
                             .create("x")
            );
        options.addOption(
                OptionBuilder.withLongOpt("rpkm-cutoff")
                             .withDescription("RPKM cutoff when using background BAM file. Genes below the cutoff will be ignored")
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
                OptionBuilder.withLongOpt("gene-selection-method")
                             .withDescription("Method to use when selecting genes for fusions: uniform|empirical|binned")
                             .hasArg()
                             .create("m")
            );
        options.addOption(
                OptionBuilder.withLongOpt("auto-correct-orientation")
                             .withDescription("Auto correct orientation of genes selected for a fusion if located on different strands")
                             .create("a")
            );
        options.addOption(
                OptionBuilder.withLongOpt("out-of-frame")
                             .withDescription("Allow fusion genes outside of reading frames")
                             .create("d")
            );
        options.addOption(
                OptionBuilder.withLongOpt("keep-exon-boundries")
                             .withDescription("Generate fusion breaks on exon boundries only")
                             .create("e")
            );
        options.addOption(
                OptionBuilder.withLongOpt("limit")
                             .withDescription("Limit fusions to specific genes/transcripts")
                             .hasArg()
                             .create("l")
            );
        options.addOption(
                OptionBuilder.withLongOpt("version")
                             .withDescription("Display version info")
                             .create("v")
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

    private void printVersionAndExit() {
        Properties properties = new Properties();
        try {
            InputStream inStream = this.getClass().getClassLoader().getResourceAsStream("version.properties");
            properties.load(inStream);
        } catch (Exception e){
            logger.warn("Failed to load version data: "+e.getMessage());
        }
        System.out.println("v"+properties.getProperty("fusim.version"));
        System.exit(0);
    }
}
