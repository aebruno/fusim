package edu.buffalo.fusim;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import edu.buffalo.fusim.gtf.GTFParseException;

public class UCSCRefFlatParser implements GeneModelParser {
    private boolean keepExonBoundries;
    private boolean cdsOnly;
    private Map<String, Boolean> limit;

    public UCSCRefFlatParser() {
        this.keepExonBoundries = false;
        this.cdsOnly = false;
        this.limit = null;
    }

    public UCSCRefFlatParser(boolean keepExonBoundries, boolean cdsOnly, Map<String,Boolean> limit) {
        this.keepExonBoundries = keepExonBoundries;
        this.cdsOnly = cdsOnly;
        this.limit = limit;
    }

    public TranscriptRecord parseLine(String line) throws GTFParseException {
        String[] fields = line.split("\t");
        TranscriptRecord record = TranscriptRecord.fromRefFlat(fields);
        
        //XXX skip the haplotypes and unassembled chroms
        if(record.getChrom().contains("_")) return null;

        if(record.getExonBases() == 0) return null;

        if(cdsOnly && record.getCdsExonBases() == 0) return null;

        if(keepExonBoundries && record.getValidExonBoundryBreaks(cdsOnly).size() == 0) {
            return null;
        }

        if(limit != null && 
          !limit.containsKey(record.getGeneId()) &&
          !limit.containsKey(record.getTranscriptId())) return null;

        return record;
    }

    //XXX test only. remove soon
    public static void main(String[] args) throws Exception {
        UCSCRefFlatParser parser = new UCSCRefFlatParser();

        FileInputStream in = new FileInputStream("../fusim-data/refFlat.txt");
        BufferedReader reader = new BufferedReader(new InputStreamReader(in));
        String line = null;
        while ((line = reader.readLine()) != null) {
            TranscriptRecord f = parser.parseLine(line);
            //System.err.println(Arrays.toString(f.generateExonBoundryBreak(true)));
            //if(f!=null && f.getCodingExons().size() == 0) System.out.println(f.getTranscriptId());
            //if(f != null && ("NR_002206".equals(f.getTranscriptId()) ||   "NR_026911".equals(f.getTranscriptId()))) {
            //    System.out.print(f);
            //    System.out.println(Arrays.toString(f.generateExonBreak(false, true)));
            //} 
        }
    }
}
