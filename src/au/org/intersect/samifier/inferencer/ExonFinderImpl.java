package au.org.intersect.samifier.inferencer;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

import au.org.intersect.samifier.domain.ExonInfo;
import au.org.intersect.samifier.domain.GenomeConstant;
import au.org.intersect.samifier.domain.MegaExonInfo;
import au.org.intersect.samifier.domain.TranscriptInfo;

public class ExonFinderImpl implements ExonFinder {

    private static Logger LOG = Logger.getLogger(ExonFinderImpl.class);

    private List<TranscriptInfo> referenceTranscripts;

    public ExonFinderImpl(List<TranscriptInfo> referenceTranscripts) {
        this.referenceTranscripts = referenceTranscripts;
    }

    public MegaExonInfo getClosestKnownExon(TranscriptInfo transcript) {
        List<ExonInfo> targetExons            = transcript.getAllExons();
        List<ExonInfo> candidateExons         = getCandidateExons(targetExons);
        List<MegaExonInfo> candidateMegaExons = convertToMegaExons(candidateExons);

        // Remove exons that do not satisfy our criteria.
        candidateMegaExons = filterCandidateExons(candidateMegaExons);

        // Iterates through the list of 'potential' exons and
        // find one with the 'best' fit. Ideally, we want an
        // exon somewhere between the 'Start' exon and the 
        // 'Stop' exon (i.e., in the middle of the transcript).
        MegaExonInfo prevExon = null;
        for (MegaExonInfo currExon : candidateMegaExons) {
            // If there are multiple 'candidate' exons, the 
            // exon with the highest score is chosen to
            // be the 'closest known'.
            prevExon = getBetterExon(prevExon, currExon);
        }

        if (prevExon == null) {
            LOG.warn("No usable exon in reference for "+ transcript.getId() 
                     + ". Must do 3/6-frame translation.");
        }
        return prevExon;
    }

    private List<ExonInfo> getCandidateExons(List<ExonInfo> targetExons) {
        List<ExonInfo> exons = new ArrayList<ExonInfo>();

        // Find exons in the reference that are
        // identical to an exon in the target
        for (TranscriptInfo t : referenceTranscripts) {
            List<ExonInfo> referenceExons = t.getAllExons();
            referenceExons.retainAll(targetExons);
            exons.addAll(referenceExons);
        }
        return exons;
    }

    private List<MegaExonInfo> convertToMegaExons(List<ExonInfo> exons) {
        // Group identical exons
        Map<ExonInfo, List<ExonInfo>> map = new HashMap<ExonInfo, List<ExonInfo>>();
        for (ExonInfo e : exons) {
            List<ExonInfo> exonList = map.get(e);
            if (exonList == null) {
                exonList = new ArrayList<ExonInfo>();
            }
            exonList.add(e);
            map.put(e, exonList);
        }

        // Merge all identical exons into a single 'Mega' exon
        List<MegaExonInfo> megaExons = new ArrayList<MegaExonInfo>();
        for (Collection<ExonInfo> eList : map.values()) {
            MegaExonInfo mExon = new MegaExonInfo((List<ExonInfo>) eList);
            megaExons.add(mExon);
        }
        return megaExons;
    }

    private List<MegaExonInfo> filterCandidateExons(List<MegaExonInfo> candidateExons) {
        Iterator<MegaExonInfo> iter = candidateExons.iterator();
        while (iter.hasNext()) {
            MegaExonInfo currExon = iter.next();
            if (!isValid(currExon)) {
                // Exons that do not satisfy our criteria will be removed
                iter.remove();
            }
        }
        return candidateExons;
    }

    private boolean isValid(MegaExonInfo exon) {
        // Must have a stop phase or else cannot use
        if (exon.getCodingFrames().isEmpty()) {
            return false;
        }

        // Exons that have more than 1 start or 
        // stop codon are considered invalid
        if (exon.getStartCodons().size() > 1
            || exon.getStopCodons().size() > 1) {
            return false;
        }

        // Exons that have a stop codon but no start 
        // codon are considered invalid beause we cannot 
        // determine which M will be the starting M
        if (exon.getStartCodons().size() == 0
            && exon.getStopCodons().size() == 1) {
            return false;
        }

        // Exons that have a start AND a stop codon that do not satify
        // the following criteria are considered invalid.
        //     1) Stop codon MUST be after the start codon.
        //     2) The nucleotide sequence MUST have at least 2 amino
        //        acids worth of coding. (Smallest valid sequence = M*)
        //     3) The nucleotide sequence MUST be a multiple of 3. (If not,
        //        there is more potential of getting the incorrect
        //        sequence in the later stages)
        if (exon.getStartCodons().size() == 1
            && exon.getStopCodons().size() == 1) {
            int startCodon = exon.getStartCodons().iterator().next();
            int stopCodon  = exon.getStopCodons().iterator().next();
            if (exon.isForward() && startCodon > stopCodon) {
                return false;

            } else if (!exon.isForward() && startCodon < stopCodon) {
                return false;
            }

            int numBases   = Math.abs(stopCodon - startCodon) + 1;
            if (numBases < (GenomeConstant.BASES_PER_CODON * 2) 
                || numBases % GenomeConstant.BASES_PER_CODON  != 0) {
                return false;
            }
        }
        return true;
    }

    public MegaExonInfo getBetterExon(MegaExonInfo prevExon, MegaExonInfo currExon) {
        if (prevExon == null) {
            return currExon;
        }

        if (currExon.getScore() > prevExon.getScore()) {
            return currExon;
        }

        if (!currExon.getStartCodons().isEmpty() 
            && prevExon.getStartCodons().isEmpty()) {
            return currExon;
        }

        if (currExon.getStopCodons().size() 
            > prevExon.getStopCodons().size()) {
            return currExon;
        }

        return prevExon;
    }

}

