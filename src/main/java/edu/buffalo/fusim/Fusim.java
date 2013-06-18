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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Random;

import net.sf.picard.sam.MergeSamFiles;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.io.IOUtils;
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

        
        int nFusions = 0;
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

        int nTriFusion = 0;
        if(cmd.hasOption("j")) {
            try {
                nTriFusion = Integer.parseInt(cmd.getOptionValue("j"));
            } catch(NumberFormatException e) {
                printHelpAndExit(options, "Number of tri-fusions (-j) must be a number");
            }
        }

        int nIntraChromFusion = 0;
        if(cmd.hasOption("y")) {
            try {
                nIntraChromFusion = Integer.parseInt(cmd.getOptionValue("y"));
            } catch(NumberFormatException e) {
                printHelpAndExit(options, "Number of intra-chromosome fusions (-y) must be a number");
            }
        }

        int nSelfFusion = 0;
        if(cmd.hasOption("s")) {
            try {
                nSelfFusion = Integer.parseInt(cmd.getOptionValue("s"));
            } catch(NumberFormatException e) {
                printHelpAndExit(options, "Number of self-fusions (-s) must be a number");
            }
        }

        int foreignInsertionLen = 0;
        if(cmd.hasOption("u")) {
            try {
                foreignInsertionLen = Integer.parseInt(cmd.getOptionValue("u"));
            } catch(NumberFormatException e) {
                printHelpAndExit(options, "Foreign insertion length (-u) must be a number");
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

        double foreignInsertionPct = 0.0;
        if(cmd.hasOption("w")) {
            try {
                foreignInsertionPct = Double.parseDouble(cmd.getOptionValue("w"));
                if(foreignInsertionPct < 0 || foreignInsertionPct > 1) throw new NumberFormatException();
            } catch(NumberFormatException e) {
                printHelpAndExit(options, "Foreign insertion percent (-w) must be 0 < x < 1");
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
            if(limits.length < 1) {
                printHelpAndExit(options, "Must provide a limit of at least one genes (ex. gene1,gene2,..): "+cmd.getOptionValue("l"));
            }
            for(int i = 0; i < limits.length; i++) {
                limit.put(limits[i], true);
            }
        }

        List<String[]> filters = new ArrayList<String[]>();
        for(String filterOption : new String[]{"1","2","3"}) {
            if(cmd.hasOption(filterOption)) {
                filters.add(cmd.getOptionValue(filterOption).split(","));
            } else {
                filters.add(null);
            }
        }
        
        File bamFile = null;
        if(cmd.hasOption("b")) {
            bamFile = new File(cmd.getOptionValue("b"));
        
            if(!bamFile.canRead()) {
                printHelpAndExit(options, "Please provide a valid BAM file");
            }
        }
        
        logger.info("========================================================================");
        logger.info("Running Fusim with the following settings:");
        logger.info("========================================================================");
        logger.info("Reference Gene Model: "+geneModelFile.getAbsolutePath());
        if(cmd.hasOption("t")) {
            logger.info("Text Output: "+("-".equals(cmd.getOptionValue("t")) ? "<stdout>" : cmd.getOptionValue("t")));
        }
        if(cmd.hasOption("f")) {
            logger.info("Fasta Output: "+("-".equals(cmd.getOptionValue("f")) ? "<stdout>" : cmd.getOptionValue("f")));
        }
        if(!cmd.hasOption("f") && !cmd.hasOption("t")) {
            logger.info("Text Output: <stdout>");
        }
        logger.info("");
        logger.info("------------------");
        logger.info("Gene Selection");
        logger.info("------------------");
        if(cmd.hasOption("b")) {
            logger.info("Mode: background reads");
            logger.info("BAM file: "+bamFile.getAbsolutePath());
            logger.info("RPKM cutoff: "+rpkmCutoff);
            logger.info("Number of threads: "+nThreads);
            logger.info("Gene selection method: "+geneSelectioMethod.toString());
        } else {
            logger.info("Mode: gene model");
        }
        logger.info("");
        logger.info("------------------");
        logger.info("Type of fusions");
        logger.info("------------------");
        logger.info("Hybrid: "+nFusions);
        logger.info("Self: "+nSelfFusion);
        logger.info("Complex: "+nTriFusion);
        logger.info("Intra-chromosome: "+nIntraChromFusion);
        logger.info("Read through: "+nReadThrough);
        logger.info("");
        logger.info("------------------");
        logger.info("Fusion options");
        logger.info("------------------");
        logger.info("CDS only: "+(cmd.hasOption("c") ? "yes" : "no"));
        logger.info("Auto-correct orientation: "+(cmd.hasOption("a") ? "yes" : "no"));
        logger.info("Allow fusions outside of ORF: "+(cmd.hasOption("d") ? "yes" : "no"));
        logger.info("Force fusion breaks on exon boundries: "+(cmd.hasOption("e") ? "yes" : "no"));
        if(cmd.hasOption("u")) {
            logger.info("Foreign insertion max length: "+foreignInsertionLen);
            logger.info("Foreign insertion percent: "+foreignInsertionPct);
        }
        logger.info("========================================================================");
        
        GeneModelParser parser = new UCSCRefFlatParser(cmd.hasOption("e"), cmd.hasOption("c"), limit);
        GeneSelector selector = null;
        FusionGenerator fg = null;
        
        if(cmd.hasOption("b")) {
            selector = new BackgroundSelector(bamFile, rpkmCutoff, nThreads);
            fg = new BackgroundGenerator();
        } else {
            selector = new StaticSelector();
            fg = new RandomGenerator();
        }

        selector.setGeneModelFile(geneModelFile);
        selector.setGeneModelParser(parser);

        fg.setGeneSelector(selector);
        fg.setGeneSelectionMethod(geneSelectioMethod);
        fg.setFilters(filters);
        
        List<FusionGene> fusions = new ArrayList<FusionGene>();
        
        if(nFusions > 0) {
            fusions.addAll(fg.generate(nFusions, 2));
        }
        
        // Generate any read through fusion genes
        if(nReadThrough > 0) {
            logger.info("Generating read through genes...");
            ReadThroughGenerator rt = new ReadThroughGenerator();
            rt.setGeneSelector(selector);
            rt.setGeneSelectionMethod(geneSelectioMethod);

            List<FusionGene> rtFusions = rt.generate(nReadThrough, 2);
            for(FusionGene g : rtFusions) {
                g.setFusionType(FusionType.READ_THROUGH);
            }
            fusions.addAll(rtFusions);
        }
        
        // Generate any tri-fusions
        if(nTriFusion > 0) {
            logger.info("Generating tri-fusion genes...");
            List<FusionGene> tfusions = fg.generate(nTriFusion, 3);
            for(FusionGene g : tfusions) {
                g.setFusionType(FusionType.TRI_FUSION);
            }
            fusions.addAll(tfusions);
        }
        
        // Generate any intra chromosome fusions
        if(nIntraChromFusion > 0) {
            logger.info("Generating intra-chromosome fusions...");
            IntraChromGenerator ig = new IntraChromGenerator();
            ig.setGeneSelector(selector);
            ig.setGeneSelectionMethod(geneSelectioMethod);

            List<FusionGene> ifusions = ig.generate(nIntraChromFusion, 2);
            for(FusionGene g : ifusions) {
                g.setFusionType(FusionType.INTRA_CHROMOSOME);
            }
            fusions.addAll(ifusions);
        }
        
        // Generate any self-fusions
        if(nSelfFusion > 0) {
            logger.info("Generating self-fusion genes...");
            List<FusionGene> sfusions = fg.generate(nSelfFusion, 1);
            for(FusionGene g : sfusions) {
                g.setFusionType(FusionType.SELF_FUSION);
            }
            fusions.addAll(sfusions);
        }

        if(fusions.size() == 0) {
            fatalError("No fusions to simulate! Check to be sure you have -j,-n,-s,-x,-y specified and your filters are correct.");    
        }
        
        if(textOutput != null) {
            textOutput.println(StringUtils.join(FusionGene.getHeader(), "\t"));
        }

        int foreignInsertionCutoff = (int)(foreignInsertionPct*fusions.size());
        
        Random rgen = new Random();
        for(int g = 0; g < fusions.size(); g++) {
            FusionGene f = fusions.get(g);

            //out.println(f);
            List<int []> breaks = new ArrayList<int []>();
            
            // First half of gene 1
            breaks.add(f.getGene(0).generateExonBreak(true, cmd.hasOption("c")));
            
            if(f.size() == 2) {
                // Second half of gene2
                breaks.add(f.getGene(1).generateExonBreak(false, cmd.hasOption("c")));
            } else if(f.size() == 3) {
                // Second half of gene2
                breaks.add(f.getGene(1).generateExonBreak(false, cmd.hasOption("c")));
                
                // Second half of gene3
                breaks.add(f.getGene(2).generateExonBreak(false, cmd.hasOption("c")));
            }

            // Keep ORF (don't allow out of frame) and allow splitting of exons
            if(!cmd.hasOption("d") && !cmd.hasOption("e")) {
                // Split last exon in half and ensure within ORF
                for(int i = 0; i < breaks.size(); i++) {
                    int[] exons = breaks.get(i);
                    int[] lastExon = f.getGene(i).getExons(cmd.hasOption("c")).get(exons[exons.length-1]);
                    int randIndex = rgen.nextInt(lastExon[1]-lastExon[0]);
                    while(randIndex % 3 != 0) {
                        randIndex--;
                    }
                    f.getGene(i).getExons(cmd.hasOption("c")).get(exons[exons.length-1])[1] -= randIndex;
                }
            } else if(cmd.hasOption("e") && !cmd.hasOption("d")) {
                breaks.clear();
                // Keep ORF (don't allow out of frame) and don't allow splitting of exons (keep exon boundries)
                // Break genes on exons boundries
                breaks.add(f.getGene(0).generateExonBoundryBreak(cmd.hasOption("c")));
                
                if(f.size() == 2) {
                    breaks.add(f.getGene(1).generateExonBoundryBreak(cmd.hasOption("c")));
                } else if(f.size() == 3) {
                    breaks.add(f.getGene(2).generateExonBoundryBreak(cmd.hasOption("c")));
                    breaks.add(f.getGene(3).generateExonBoundryBreak(cmd.hasOption("c")));
                }
            }

            // Set options for output
            if(cmd.hasOption("a")) {
                f.addOption(FusionOption.AUTO_CORRECT_ORIENTATION);
            }
            if(cmd.hasOption("c")) {
                f.addOption(FusionOption.CDS_ONLY);
            }
            if(cmd.hasOption("d")) {
                f.addOption(FusionOption.OUT_OF_FRAME);
            } else {
                f.addOption(FusionOption.SYMMETRICAL_EXONS);
            }
            if(cmd.hasOption("e")) {
                f.addOption(FusionOption.KEEP_EXON_BOUNDRY);
            }
            
            if(textOutput != null) {
                textOutput.print(f.outputText(breaks, cmd.hasOption("c")));
            }
            
            if(fastaOutput != null) {
                if(foreignInsertionLen > 0 && foreignInsertionCutoff > 0 && g <= foreignInsertionCutoff) {
                    f.addOption(FusionOption.FOREIGN_INSERTION);
                    fastaOutput.println(f.outputFasta(breaks, referenceFile, cmd.hasOption("c"), cmd.hasOption("a"), foreignInsertionLen));
                } else {
                    fastaOutput.println(f.outputFasta(breaks, referenceFile, cmd.hasOption("c"), cmd.hasOption("a")));
                }
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
                         .withDescription("Path to gene model file in refFlat format")
                         .hasArg()
                         .create("g")
        );
        options.addOption(
                OptionBuilder.withLongOpt("reference")
                             .withDescription("Path to indexed reference genome fasta file (.fai)")
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
                             .withDescription("Number of fusions to generate using two randomly selected genes")
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
                             .withDescription("Number of read through fusions")
                             .hasArg()
                             .create("x")
            );
        options.addOption(
                OptionBuilder.withLongOpt("tri-fusion")
                             .withDescription("Number of fusions with three genes")
                             .hasArg()
                             .create("j")
            );
        options.addOption(
                OptionBuilder.withLongOpt("self-fusion")
                             .withDescription("Number of self-fusions (fusions with single gene)")
                             .hasArg()
                             .create("s")
            );
        options.addOption(
                OptionBuilder.withLongOpt("intra-chrom")
                             .withDescription("Number of intra-chromosome fusions (fusions within single chrom)")
                             .hasArg()
                             .create("y")
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
                             .withDescription("Method to use when selecting genes for fusions: uniform|binned|empirical|empirical-sturges")
                             .hasArg()
                             .create("m")
            );
        options.addOption(
                OptionBuilder.withLongOpt("auto-correct-orientation")
                             .withDescription("Auto correct orientation of fusion sequence if genes are located on different strands")
                             .create("a")
            );
        options.addOption(
                OptionBuilder.withLongOpt("out-of-frame")
                             .withDescription("Allow fusions outside of reading frames. By default the reading frame is preserved")
                             .create("d")
            );
        options.addOption(
                OptionBuilder.withLongOpt("keep-exon-boundary")
                             .withDescription("Generate fusion breaks on exon boundaries only")
                             .create("e")
            );
        options.addOption(
                OptionBuilder.withLongOpt("limit")
                             .withDescription("Limit all fusions to specific geneId, transcriptId, or chrom")
                             .hasArg()
                             .create("l")
            );
        options.addOption(
                OptionBuilder.withLongOpt("gene1")
                             .withDescription("Filter for gene1")
                             .hasArg()
                             .create("1")
            );
        options.addOption(
                OptionBuilder.withLongOpt("gene2")
                             .withDescription("Filter for gene2")
                             .hasArg()
                             .create("2")
            );
        options.addOption(
                OptionBuilder.withLongOpt("gene3")
                             .withDescription("Filter for gene3")
                             .hasArg()
                             .create("3")
            );
        options.addOption(
                OptionBuilder.withLongOpt("foreign-insertion-length")
                             .withDescription("Maxium length of randomly generated sequence to insert between fusion breakpoints")
                             .hasArg()
                             .create("u")
            );
        options.addOption(
                OptionBuilder.withLongOpt("foreign-insertion-perecent")
                             .withDescription("Percent of fusions to insert foreign sequence between fusion breakpoints")
                             .hasArg()
                             .create("w")
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
        // Print CLI options
        //HelpFormatter formatter = new HelpFormatter();
        //formatter.printHelp("fusim", options);
        try {
            System.out.print(
                IOUtils.toString(this.getClass().getClassLoader().getResourceAsStream("fusim.options"))
            );
        } catch (Exception e){
            logger.fatal("Failed to load options help file: "+e.getMessage());
        }
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
            logger.fatal("Failed to load version data: "+e.getMessage());
        }
        System.out.println("v"+properties.getProperty("fusim.version"));
        System.exit(0);
    }
}
