package au.org.intersect.samifier.domain;

import java.util.Iterator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;
import org.apache.log4j.Logger;

public class ProteinOutputter implements Outputter {

    public static final int FASTA_LINE_LENGTH = 60;
    private static Logger LOG = Logger.getLogger(ProteinOutputter.class);

    private String fastaHeader;
    private ProteinLocation proteinLocation;
    private StringBuilder genomeString;
    private CodonTranslationTable translationTable;

    public ProteinOutputter(ProteinLocation proteinLocation,
        String databaseName, StringBuilder genomeString,
        CodonTranslationTable translationTable) {
        this.fastaHeader      = ">gn1|" + databaseName + "|" + proteinLocation.getName();
        this.proteinLocation  = proteinLocation;
        this.genomeString     = genomeString;
        this.translationTable = translationTable;
    }

    @Override
    public String getOutput() throws OutputException {
        String lineFeed = System.getProperty("line.separator");

        StringBuilder buffer = new StringBuilder();
        buffer.append(fastaHeader);
        buffer.append(lineFeed);

        int startIndex = proteinLocation.getStartIndex() - 1;
        int stopIndex = startIndex + proteinLocation.getLength();
        String sequence = genomeString.substring(startIndex, stopIndex);

        String aminoAcidSequence = null;
        if (proteinLocation.getDirection().equals(GenomeConstant.REVERSE_FLAG)) {
            StringBuilder invertedReversedSequence = new StringBuilder(
                    invertNucleotideSequence(sequence.toString())).reverse();
            aminoAcidSequence = translationTable
                    .proteinToAminoAcidSequence(invertedReversedSequence
                            .toString());
        } else {
            aminoAcidSequence = translationTable
                    .proteinToAminoAcidSequence(sequence.toString());
        }
        int sequenceLength = aminoAcidSequence.length();
        int wholeParts = sequenceLength / FASTA_LINE_LENGTH;
        int sequenceCursor = 0;
        for (int i = 0; i < wholeParts; i++) {
            buffer.append(aminoAcidSequence.substring(sequenceCursor,
                    sequenceCursor + FASTA_LINE_LENGTH));
            buffer.append(lineFeed);
            sequenceCursor += FASTA_LINE_LENGTH;
        }
        if (sequenceCursor < sequenceLength) {
            buffer.append(aminoAcidSequence.substring(sequenceCursor,
                    sequenceLength));
            buffer.append(lineFeed);
        }
        return buffer.toString();
    }

    public static String invertNucleotideSequence(String sequence) {
        return StringUtils.replaceChars(sequence, "ACGT", "TGCA");
    }

    // For generating the protein sequence of a transcript isoform
    // Constructor
    private TranscriptInfo transcript;

    public ProteinOutputter(TranscriptInfo transcript,
            String databaseName, StringBuilder genomeString,
            CodonTranslationTable translationTable) {

        this.fastaHeader = ">gn1|" + databaseName + "|" + transcript.getId();
        this.transcript = transcript;
        this.genomeString = genomeString;
        this.translationTable = translationTable;
    }

    // Functions
    public String getTranscriptOutput()
            throws OutputException {

        try {
            String startEId    = transcript.getStartCodon().getExonId();
            ExonInfo startExon = transcript.getExon(startEId);
            int startExonIdx   = transcript.getAllExons().indexOf(startExon);

            String stopEId    = transcript.getStopCodon().getExonId();
            ExonInfo stopExon = transcript.getExon(stopEId);
            int stopExonIdx   = transcript.getAllExons().indexOf(stopExon) + 1;
            
            List<String> exonSequences = transcript.getExonSequences(genomeString);
            exonSequences = exonSequences.subList(startExonIdx, stopExonIdx);
            String nucleotideSequence = String.join("", exonSequences);
            String aminoAcidSequence  = translationTable.proteinToAminoAcidSequence(nucleotideSequence);

            if (!aminoAcidSequence.startsWith("M")) {
                LOG.warn("Amino acid sequence of Transcript " + transcript.getId()
                    + " does not start with M.");
            }

            if (!aminoAcidSequence.endsWith("*")) {
                LOG.warn("Amino acid sequence of Transcript " + transcript.getId()
                    + " does not end with \\*.");
            }

            String fastaSequence = getFastaSequence(aminoAcidSequence);
            return fastaSequence;

        } catch (UnknownCodonException e) {
            String line1 = "Failed translation of Transcript " + transcript.getId() + ".";
            String line2 = "> " + e + "\n";
            String message = line1 + line2;
            throw new OutputException(message);
        }
    }

    private String getFastaSequence(String aminoAcidSequence) {
        String lineFeed      = System.getProperty("line.separator");
        int sequenceLength   = aminoAcidSequence.length();
        int wholeParts       = sequenceLength / FASTA_LINE_LENGTH;
        int sequenceCursor   = 0;

        StringBuilder buffer = new StringBuilder();
        buffer.append(fastaHeader);
        buffer.append(lineFeed);

        for (int i = 0; i < wholeParts; i++) {
            buffer.append(aminoAcidSequence.substring(sequenceCursor,
                    sequenceCursor + FASTA_LINE_LENGTH));
            buffer.append(lineFeed);
            sequenceCursor += FASTA_LINE_LENGTH;
        }
        if (sequenceCursor < sequenceLength) {
            buffer.append(aminoAcidSequence.substring(sequenceCursor,
                    sequenceLength));
            buffer.append(lineFeed);
        }

        return buffer.toString();
    }

}
