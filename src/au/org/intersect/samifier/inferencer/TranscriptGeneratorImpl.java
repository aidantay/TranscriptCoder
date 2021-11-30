package au.org.intersect.samifier.inferencer;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

import au.org.intersect.samifier.parser.FastaParser;
import au.org.intersect.samifier.domain.CodonInfo;
import au.org.intersect.samifier.domain.CodonTranslationTable;
import au.org.intersect.samifier.domain.ExonInfo;
import au.org.intersect.samifier.domain.GenomeConstant;
import au.org.intersect.samifier.domain.MegaExonInfo;
import au.org.intersect.samifier.domain.ProteinOutputter;
import au.org.intersect.samifier.domain.TranscriptomeConstant;
import au.org.intersect.samifier.domain.TranscriptInfo;
import au.org.intersect.samifier.domain.Transcriptome;
import au.org.intersect.samifier.domain.UnknownCodonException;

public class TranscriptGeneratorImpl implements TranscriptGenerator {

    private static Logger LOG = Logger.getLogger(TranscriptGeneratorImpl.class);

    private StringBuilder genomeString;
    private CodonTranslationTable translationTable;

    public TranscriptGeneratorImpl(StringBuilder genomeString, 
        CodonTranslationTable translationTable) {

        this.genomeString     = genomeString;
        this.translationTable = translationTable;
    }

    public List<TranscriptInfo> inferTranscript(TranscriptInfo prevTranscript, 
        MegaExonInfo closestKnownExon) 
            throws TranscriptGeneratorException {

        List<TranscriptInfo> inferredTranscripts = new ArrayList<TranscriptInfo>();

        // Special case where the closest known exon contains
        // both start and stop codons, and only has exactly one coding frame
        if (!closestKnownExon.getStartCodons().isEmpty()
            && !closestKnownExon.getStopCodons().isEmpty()
            && closestKnownExon.getCodingFrames().size() == 1) {

            TranscriptInfo currTranscript = getTranscripts(prevTranscript, closestKnownExon);
            inferredTranscripts.add(currTranscript);

        } else {
            // Create a new transcript for every coding frame the closest known
            // exon has (Except when there is no start codon in the transcript)
            for (Integer codingFrame : closestKnownExon.getCodingFrames()) {
                TranscriptInfo currTranscript = getTranscripts(prevTranscript, closestKnownExon, codingFrame);
                if (currTranscript == null) {
                    continue;
                }

                inferredTranscripts.add(currTranscript);
            }
        }

        if (inferredTranscripts.size() == 0) {
            LOG.warn("No appropriate coding frame for Transcript " + prevTranscript.getId()
                + ". Must do 3/6-frame translation.");

        } else {
            // Update the ID of inferred transcripts 
            int counter = 1;
            for (TranscriptInfo t : inferredTranscripts) {
                String newId = t.getId() + "_" + String.valueOf(counter);
                t.setId(newId);
                counter++;
            }
        }

        return inferredTranscripts;
    }

    private TranscriptInfo getTranscripts(TranscriptInfo prevTranscript,
        MegaExonInfo closestKnownExon) {

        // Modify our new Transcript based
        // on the information we've gathered
        TranscriptInfo currTranscript = new TranscriptInfo(prevTranscript);

        ExonInfo middleExon = prevTranscript.getExon((ExonInfo) closestKnownExon);
        int prevExonIdx = prevTranscript.getAllExons().indexOf(middleExon) - 1;
        int nextExonIdx = prevTranscript.getAllExons().indexOf(middleExon) + 1;

        // Infer the start and stop codons based on the 
        // position of the start and stop codons
        int startCodon = closestKnownExon.getStartCodons().iterator().next();
        int stopCodon  = closestKnownExon.getStopCodons().iterator().next();

        // Count the number of bases between the start and stop codons
        // and the ends of the exon
        int startUtrLength = startCodon - middleExon.getStart();
        int stopUtrLength  = middleExon.getStop() - stopCodon;
        if (!currTranscript.isForward()) {
            startUtrLength = middleExon.getStop() - startCodon;
            stopUtrLength  = stopCodon - middleExon.getStart();
        }

        // Account for exons (if any) before and after the middle exon
        if (prevExonIdx > -1) {
            ExonInfo tempMiddleExon = currTranscript.getAllExons().get(prevExonIdx);
            String startNucSequence = getStartNucSequence(currTranscript, tempMiddleExon);
            startUtrLength = startUtrLength + startNucSequence.length();
        }

        if (nextExonIdx < prevTranscript.getAllExons().size()) {
            String stopNucSequence = getStopNucSequence(prevTranscript, nextExonIdx);
            stopUtrLength = stopUtrLength + stopNucSequence.length();
        }

        // Modify our new Transcript based on the information we've gathered
        modifyStartTranscript(currTranscript, startUtrLength);
        modifyStopTranscript(currTranscript, stopUtrLength);
        return currTranscript;
    }

    private TranscriptInfo getTranscripts(TranscriptInfo prevTranscript,
        MegaExonInfo closestKnownExon, int codingFrame) 
            throws TranscriptGeneratorException {

        try {
            // Create a new transcript that we can modifier
            TranscriptInfo currTranscript = new TranscriptInfo(prevTranscript);
            ExonInfo middleExon = prevTranscript.getExon((ExonInfo) closestKnownExon);

            // Infer the start by translating the nucleotide sequence
            // and calculating the UTR
            String startNucSequence = getStartNucSequence(currTranscript, middleExon); 
            int startUtrLength      = getStartUtrLength(startNucSequence, codingFrame);
            // If we cannot find a start codon, then we
            // assume the translation frame is invalid
            // and try another coding frame.
            if (startUtrLength == -1) {
                return null;
            }

            // Modify our new Transcript based
            // on the information we've gathered
            modifyStartTranscript(currTranscript, startUtrLength);

            // Infer the stop by translating the nucleotide sequence
            // and calculating the UTR
            String stopNucSequence = getStopNucSequence(currTranscript); 
            int stopUtrLength      = getStopUtrLength(stopNucSequence);

            // Modify our new Transcript based
            // on the information we've gathered
            modifyStopTranscript(currTranscript, stopUtrLength);
            return currTranscript;

        } catch (TranscriptGeneratorException e) {
            String message = "Could not find a " + e + " for Transcript " + prevTranscript.getId()
                + ". Failed translation.";
            throw new TranscriptGeneratorException(message);
        }
    }

    private String getStartNucSequence(TranscriptInfo transcript, 
        ExonInfo middleExon) {

        int startExonIdx  = 0;
        int middleExonIdx = transcript.getAllExons().indexOf(middleExon);

        List<String> exonSequences      = transcript.getExonSequences(genomeString);
        List<String> startExonSequences = exonSequences.subList(startExonIdx, middleExonIdx + 1);
        String startNucSequence = String.join("", startExonSequences);
        return startNucSequence;
    }

    private int getStartUtrLength(String startNucSequence, int codingFrame) 
            throws TranscriptGeneratorException {

        try {
            /**********
             * Before we do anything, we need to fix
             * the nucleotide sequence for translation
             **********/
            // Remove bases from the end (according to the exon coding frame)
            int stopBases    = startNucSequence.length() - codingFrame;
            startNucSequence = startNucSequence.substring(0, stopBases);

            // Remove bases from the beginning (to fit translation)
            int startBases   = startNucSequence.length() % GenomeConstant.BASES_PER_CODON;
            startNucSequence = startNucSequence.substring(startBases);

            // Translate the nucleotide sequence
            String startAaSequence    = translationTable.proteinToAminoAcidSequence(startNucSequence);
            String newStartAaSequence = getStartAaSequence(startAaSequence);
            if (newStartAaSequence == null) {
                // If there aren't any 'M', then there aren't any start sites.
                // We therefore assume that there isn't any UTR.
                return -1;
            }

            // Calculate the UTR based on differences in the nucleotide
            // and amino acid sequences
            int seqDiff    = startAaSequence.length() - newStartAaSequence.length();
            int startUtrLength = startBases + (seqDiff * 3);
            return startUtrLength;

        } catch (UnknownCodonException e) {
            String message = TranscriptomeConstant.START_CODON;
            throw new TranscriptGeneratorException(message);
        }
    }

    private String getStartAaSequence(String prevAaSequence) {
        String currAaSequence = prevAaSequence;

        // Find the index of the last '*' by working backward
        // i.e., Reverse the sequence and inverse the idx 
        currAaSequence = new StringBuilder(currAaSequence).reverse().toString();
        int idx = currAaSequence.indexOf('*');
        if (idx == -1) {
            // If there aren't any '*', then we're looking
            // at a sequence containing ONLY amino acids.
            idx = 0;
        } else {
            idx = prevAaSequence.length() - idx - 1;
        }

        // Find the index of the first 'M' by working forward
        // Sequences with more than one 'M' cannot be distinguished, hence
        // we take the longest sequence 
        currAaSequence = prevAaSequence.substring(idx);
        idx = currAaSequence.indexOf('M');
        if (idx == -1) {
            // If there aren't any 'M', then there aren't any start sites.
            // We therefore assume that there isn't any UTR.
            return null;
        }

        currAaSequence = currAaSequence.substring(idx);
        return currAaSequence;
    }

    private void modifyStartTranscript(TranscriptInfo transcript, 
        int prevUtrLength) {

        CodonInfo startCodon = getStartCodon(transcript, prevUtrLength);
        transcript.addCodon(startCodon);

        // Find the exons containing the start and stop 
        // codons and split the exon (if necessary)
        if (transcript.isForward()) {
            ExonInfo startExon = transcript.getExon(startCodon.getExonId());
            if (startExon.getStart() != startCodon.getStart()) {
                float eId    = Float.parseFloat(startExon.getId());
                int utrStop  = startCodon.getStart() - 1;
                String utrId = String.valueOf(eId - 1 + 0.1);

                ExonInfo utrExon = new ExonInfo(startExon);
                utrExon.setStop(utrStop);
                utrExon.setId(utrId);
                transcript.addExon(utrExon);
            }
            startExon.setStart(startCodon.getStart());
            startExon.setStartCodon(startCodon.getStart());

        } else {
            ExonInfo startExon = transcript.getExon(startCodon.getExonId());
            if (startExon.getStop() != startCodon.getStop()) {
                float eId    = Float.parseFloat(startExon.getId());
                int utrStart = startCodon.getStop() + 1;
                String utrId = String.valueOf(eId - 1 +  0.1);

                ExonInfo utrExon = new ExonInfo(startExon);
                utrExon.setStart(utrStart);
                utrExon.setId(utrId);
                transcript.addExon(utrExon);
            }
            startExon.setStop(startCodon.getStop());
            startExon.setStartCodon(startCodon.getStop());
        }
    }

    private CodonInfo getStartCodon(TranscriptInfo transcript, 
        int prevUtrLength) {

        CodonInfo startCodon = null;
        for (ExonInfo e : transcript.getAllExons()) {
            int exonLength    = e.getStop() - e.getStart() + 1;
            int currUtrLength = prevUtrLength - exonLength;

            if (currUtrLength < 0) {
                if (e.isForward()) {
                    startCodon = createLeftCodon(e, prevUtrLength,
                        TranscriptomeConstant.START_CODON);

                } else {
                    startCodon = createRightCodon(e, prevUtrLength,
                        TranscriptomeConstant.START_CODON);
                }
                break;

            } else {
                prevUtrLength = currUtrLength;
            }
        }
        return startCodon;
    }

    private String getStopNucSequence(TranscriptInfo transcript) {
        String startEId    = transcript.getStartCodon().getExonId();
        ExonInfo startExon = transcript.getExon(startEId);

        int startExonIdx  = transcript.getAllExons().indexOf(startExon);
        int stopExonIdx   = transcript.getAllExons().size();

        List<String> exonSequences     = transcript.getExonSequences(genomeString);
        List<String> stopExonSequences = exonSequences.subList(startExonIdx, stopExonIdx);
        String stopNucSequence = String.join("", stopExonSequences);
        return stopNucSequence;
    }

    private String getStopNucSequence(TranscriptInfo transcript,
        int startExonIdx) {

        int stopExonIdx = transcript.getAllExons().size();

        List<String> exonSequences     = transcript.getExonSequences(genomeString);
        List<String> stopExonSequences = exonSequences.subList(startExonIdx, stopExonIdx);
        String stopNucSequence = String.join("", stopExonSequences);
        return stopNucSequence;
    }

    private int getStopUtrLength(String stopNucSequence) 
            throws TranscriptGeneratorException {

        try {
            /**********
             * Before we do anything, we need to fix
             * the nucleotide sequence for translation
             **********/
            // Remove bases from the end (according to the exon coding frame)
            int endBases    = (stopNucSequence.length() % GenomeConstant.BASES_PER_CODON);
            stopNucSequence = stopNucSequence.substring(0, stopNucSequence.length() - endBases);

            // Translate the nucleotide sequence
            String stopAaSequence    = translationTable.proteinToAminoAcidSequence(stopNucSequence);
            String newStopAaSequence = getStopAaSequence(stopAaSequence);

            int seqDiff   = stopAaSequence.length() - newStopAaSequence.length();
            int stopUtrLength = endBases + (seqDiff * 3);
            return stopUtrLength;

        } catch (UnknownCodonException e) {
            String message = TranscriptomeConstant.STOP_CODON;
            throw new TranscriptGeneratorException(message);
        }
    }

    private String getStopAaSequence(String prevAaSequence) {
        String currAaSequence = prevAaSequence;

        // Theres a potential to find a '*' that we don't actually want
        // towards the beginning of the sequence. Not exactly sure how
        // we can handle these cases (without having to shuffle the 
        // code a bit).

        // Find the index of the first '*'
        // We stop at the first stop codon we encounter
        int idx = currAaSequence.indexOf('*');
        if (idx == -1) {
            // If there aren't any '*', then there aren't any stop sites.
            // We therefore assume that there isn't any UTR.
            idx = currAaSequence.length() - 1;
        }
        currAaSequence = currAaSequence.substring(0, idx + 1);
        return currAaSequence;
    }

    private void modifyStopTranscript(TranscriptInfo transcript,
        int prevUtrLength) {

        CodonInfo stopCodon = getStopCodon(transcript, prevUtrLength);
        transcript.addCodon(stopCodon);

        // Find the exons containing the start and stop 
        // codons and split the exon (if necessary)
        if (transcript.isForward()) {
            ExonInfo stopExon = transcript.getExon(stopCodon.getExonId());
            if (stopExon.getStop() != stopCodon.getStop()) {
                float eId    = Float.parseFloat(stopExon.getId());
                int utrStart = stopCodon.getStop() + 1;
                String utrId = String.valueOf(eId + 0.1);

                ExonInfo utrExon = new ExonInfo(stopExon);
                utrExon.setStart(utrStart);
                utrExon.setId(utrId);
                transcript.addExon(utrExon);
            }
            stopExon.setStop(stopCodon.getStop());
            stopExon.setStopCodon(stopCodon.getStop());

        } else {
            ExonInfo stopExon = transcript.getExon(stopCodon.getExonId());
            if (stopExon.getStart() != stopCodon.getStart()) {
                float eId    = Float.parseFloat(stopExon.getId());
                int utrStop  = stopCodon.getStart() - 1;
                String utrId = String.valueOf(eId + 0.1);

                ExonInfo utrExon = new ExonInfo(stopExon);
                utrExon.setStop(utrStop);
                utrExon.setId(utrId);
                transcript.addExon(utrExon);
            }
            stopExon.setStart(stopCodon.getStart());
            stopExon.setStopCodon(stopCodon.getStart());
        }
    }

    private CodonInfo getStopCodon(TranscriptInfo transcript,
        int prevUtrLength) {

        CodonInfo stopCodon  = null;
        List<ExonInfo> exons = transcript.getAllExons();
        Collections.reverse(exons);

        for (ExonInfo e : exons) {
            int exonLength    = e.getStop() - e.getStart() + 1;
            int currUtrLength = prevUtrLength - exonLength;

            if (currUtrLength < 0) {
                if (e.isForward()) {
                    stopCodon = createRightCodon(e, prevUtrLength,
                        TranscriptomeConstant.STOP_CODON);

                } else {
                    stopCodon = createLeftCodon(e, prevUtrLength,
                        TranscriptomeConstant.STOP_CODON);
                }
                break;

            } else {
                prevUtrLength = currUtrLength;
            }
        }
        return stopCodon;
    }

    private CodonInfo createLeftCodon(ExonInfo exon, 
        int prevUtrLength, String codonType) {

        int codonStartPos = exon.getStart() + prevUtrLength;
        int codonStopPos  = (codonStartPos + GenomeConstant.BASES_PER_CODON) - 1;
        if (codonStopPos > exon.getStop()) {
            codonStopPos = exon.getStop();
        }
        CodonInfo leftCodon = new CodonInfo(exon.getTranscriptId(), 
            exon.getId(), codonType, exon.getChromosome(), 
            codonStartPos, codonStopPos, exon.getDirection());
        return leftCodon;
    }

    private CodonInfo createRightCodon(ExonInfo exon, 
        int prevUtrLength, String codonType) {

        int codonStopPos  = exon.getStop() - prevUtrLength;
        int codonStartPos = (codonStopPos - GenomeConstant.BASES_PER_CODON) + 1;
        if (codonStartPos < exon.getStart()) {
            codonStartPos = exon.getStart();
        }
        CodonInfo rightCodon = new CodonInfo(exon.getTranscriptId(),
            exon.getId(), codonType, exon.getChromosome(),
            codonStartPos, codonStopPos, exon.getDirection());
        return rightCodon;
    }

}
