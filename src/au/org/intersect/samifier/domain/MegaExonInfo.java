package au.org.intersect.samifier.domain;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.builder.EqualsBuilder;
import org.apache.commons.lang3.builder.HashCodeBuilder;

public class MegaExonInfo extends ExonInfo {

    private Set<Integer> startCodons;               // List of all startCodon positions [For Transcriptome] 
    private Set<Integer> stopCodons;                // List of all stopCodon positions [For Transcriptome]
    private Set<Integer> codingFrames;              // List of all phases for this exon [For Transcriptome]

    /**********
     * Mega exons are ExonInfo objects that summarise
     * info in multiple ExonInfo objects
     **********/
    public MegaExonInfo(List<ExonInfo> exonList) {
        super(exonList.get(0));

        this.startCodons  = new HashSet<Integer>();
        this.stopCodons   = new HashSet<Integer>();
        this.codingFrames = new HashSet<Integer>();
        init(exonList);
    }

    private void init(List<ExonInfo> exonList) {
        for (ExonInfo e : exonList) {
            if (e.hasScore()) {
                if (!hasScore()) {
                    setScore(e.getScore());
                } else {
                    setScore(getScore() + e.getScore());
                }
            }

            if (e.hasCodingFrame()) {
                this.codingFrames.add(e.getCodingFrame());
            }

            if (e.containsStartCodon()) {
                this.startCodons.add(e.getStartCodon());
            }

            if (e.containsStopCodon()) {
                this.stopCodons.add(e.getStopCodon());
            }
        }
    }

    // Get Functions
    public Set<Integer> getCodingFrames() {
        return codingFrames;
    }

    public Set<Integer> getStartCodons() {
        return startCodons;
    }

    public Set<Integer> getStopCodons() {
        return stopCodons;
    }

    public String toString() {
        StringBuffer out = new StringBuffer();
        out.append(getTranscriptId() + "\t" + getId() + "\t" + getChromosome());
        out.append("\t");
        out.append(getDirectionStr() + "\t" + getStart() + "\t" + getStop());

        if (hasScore()) {
            out.append("\t");
            out.append("Score: " + getScore());
        }

        if (!codingFrames.isEmpty()) {
            out.append("\t");
            out.append("CodingFrames: " + codingFrames);
        }

        if (!startCodons.isEmpty()) {
            out.append("\t");
            out.append("Start codons: " + startCodons);
        }

        if (!stopCodons.isEmpty()) {
            out.append("\t");
            out.append("Stop codons: " + stopCodons);
        }

        return out.toString();
    }

}
