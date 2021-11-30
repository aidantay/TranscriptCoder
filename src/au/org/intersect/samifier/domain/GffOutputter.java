package au.org.intersect.samifier.domain;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

public class GffOutputter implements Outputter {
    private String genomeFileName;
    private Integer start;
    private Integer end;
    private String glimmerScore;
    private String directionFlag;
    private String glimmerName;
    private Set<String> virtualProteinNames;
    private String origin;

    public GffOutputter(ProteinLocation location,
            String genomeFileNameWithExtension) {
        this.genomeFileName = location.getChromosome(); //genomeFileNameNoExtension(genomeFileNameWithExtension);
        this.start = location.getStartIndex();
        this.end = location.getStop();
        if (location.getConfidenceScore() != null) {
            this.glimmerScore = location.getConfidenceScore().toString();
        } else {
            this.glimmerScore = "0";
        }
        this.directionFlag = location.getDirection();
        this.glimmerName = location.getName();
        this.virtualProteinNames = location.getVirtualProteinNames();
        this.origin = location.getOrigin();
    }

    @Override
    public String getOutput() {
        StringBuilder output = new StringBuilder();

        output.append(getOutputLine("gene", false));
        output.append(System.getProperty("line.separator"));

        output.append(getOutputLine("CDS", true));
        output.append(System.getProperty("line.separator"));

        return output.toString();

    }

    private String getOutputLine(String type, boolean addParent) {
        StringBuilder output = new StringBuilder();
        output.append(genomeFileName);
        if (origin != null && origin.length() >0){
            column(output, origin);
        } else {
            column(output, "Glimmer");
        }
        column(output, type);
        column(output, Integer.toString(start));
        column(output, Integer.toString(end));
        column(output, glimmerScore);
        column(output, directionFlag);
        //column(output, frame);
        column(output, "0"); // this is phase - we  currently don't use that
        ArrayList<String> attributes = new ArrayList<String>();
        attributes.add("Name=" + glimmerName);
        if (addParent) {
            attributes.add("Parent=" + glimmerName);
        } else {
            attributes.add("ID=" + glimmerName);
        }
        if (virtualProteinNames.size() > 0) {
            String virtualProteins = "";
            for (String s:virtualProteinNames) {
                virtualProteins += (virtualProteins == "" ? "" : ",") + s;
            }
            attributes.add("Virtual_protein=" + virtualProteins);
        }
        column(output, attributes);
        return output.toString();
    }

    private void column(StringBuilder buff, String field) {
        buff.append('\t');
        buff.append(field);
    }

    private void column(StringBuilder buff, List<String> attributes) {
        buff.append('\t');
        for (String attribute : attributes) {
            buff.append(attribute + ";");
        }
    }

    // For generating the GFF lines of a transcript isoform
    // Constructor
    private TranscriptInfo transcript;

    public GffOutputter(TranscriptInfo transcript) {
        this.transcript     = transcript;
        this.genomeFileName = transcript.getChromosome();
        this.start          = transcript.getStart();
        this.end            = transcript.getStop();
        this.glimmerScore   = "0";
        this.directionFlag  = transcript.getDirectionStr();
        this.glimmerName    = transcript.getId();
        this.origin         = "TranscriptCoder";
    }

    // Functions
    public String getTranscriptOutput()
            throws OutputException {

        StringBuilder output = new StringBuilder();
        output.append(getTranscriptOutputLine("gene", false));
        output.append(System.getProperty("line.separator"));

        /**********
         * Iterate through the exons and append the info
         **********/
        Iterator<ExonInfo> iter = transcript.getAllExons().iterator();
        CodonInfo startCodon    = transcript.getStartCodon();
        CodonInfo stopCodon     = transcript.getStopCodon();

        // Only need to update transcripts that have
        // Start and Stop codons
        if (startCodon != null && stopCodon != null) {
            ExonInfo prevExon = null;
            while (iter.hasNext()) {
                ExonInfo currExon = iter.next();
                float eId      = Float.parseFloat(currExon.getId());
                float startEId = Float.parseFloat(startCodon.getExonId()); 
                float stopEId  = Float.parseFloat(stopCodon.getExonId()); 

                if (prevExon != null) {
                    output.append(getExonOutputLine("intron", currExon, prevExon));
                    output.append(System.getProperty("line.separator"));
                }

                // Exon is a Start/Stop boundary
                if (eId == startEId || eId == stopEId) {
                    output.append(getExonOutputLine("CDS", currExon, null));
                    output.append(System.getProperty("line.separator"));

                // Exon is completely within the coing region (i.e., coding) 
                } else if (eId > startEId && eId < stopEId) {
                    output.append(getExonOutputLine("CDS", currExon, null));
                    output.append(System.getProperty("line.separator"));

                // Exon is completely outside the coding region (i.e., non-coding)
                } else if (eId < startEId) {
                    output.append(getExonOutputLine("five_prime_UTR_intron", currExon, null));
                    output.append(System.getProperty("line.separator"));

                } else if (eId > stopEId) {
                    output.append(getExonOutputLine("three_prime_UTR_intron", currExon, null));
                    output.append(System.getProperty("line.separator"));

                } else {
                    String errorMessage = "Exon " + eId + " in transcript " + transcript.getId() + "has unexpected type.";
                    throw new OutputException(errorMessage);
                }
                prevExon = currExon;
            }
        }
        return output.toString();
    }

    private String getExonOutputLine(String type, ExonInfo currExon, 
        ExonInfo prevExon) {

        StringBuilder output = new StringBuilder();
        output.append(genomeFileName);
        column(output, origin);
        column(output, type);

        if (prevExon == null) {
            column(output, Integer.toString(currExon.getStart()));
            column(output, Integer.toString(currExon.getStop()));
            column(output, glimmerScore);
            column(output, currExon.getDirectionStr());
            column(output, "0"); // this is phase - we  currently don't use that
            ArrayList<String> attributes = new ArrayList<String>();
            attributes.add("Name=" + glimmerName);
            attributes.add("Parent=" + glimmerName);
            column(output, attributes);
        } else {
            int start = prevExon.getStop() + 1;
            int stop  = currExon.getStart() - 1;
            if (!currExon.isForward()) {
                start = currExon.getStop() + 1;
                stop  = prevExon.getStart() - 1;
            }

            column(output, Integer.toString(start));
            column(output, Integer.toString(stop));
            column(output, glimmerScore);
            column(output, currExon.getDirectionStr());
            column(output, "0"); // this is phase - we  currently don't use that
            ArrayList<String> attributes = new ArrayList<String>();
            attributes.add("Name=" + glimmerName);
            attributes.add("Parent=" + glimmerName);
            column(output, attributes);
        }

        return output.toString();
    }

    private String getTranscriptOutputLine(String type, boolean addParent) {
        StringBuilder output = new StringBuilder();

        output.append(genomeFileName);
        column(output, origin);
        column(output, type);
        column(output, Integer.toString(start));
        column(output, Integer.toString(end));
        column(output, glimmerScore);
        column(output, directionFlag);
        column(output, "0"); // this is phase - we  currently don't use that
        ArrayList<String> attributes = new ArrayList<String>();
        attributes.add("Name=" + glimmerName);
        attributes.add("ID=" + glimmerName);
        column(output, attributes);

        return output.toString();
    }

}
