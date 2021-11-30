package au.org.intersect.samifier.domain;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.builder.EqualsBuilder;
import org.apache.commons.lang3.builder.HashCodeBuilder;

public class ExonInfo {

    private String transcriptId;
    private String exonId;
    private String chromosome;
    private int start;              // Start position of the exon (WRT Genome)
    private int stop;               // End position of the exon   (WRT Genome)
    private int direction;

    private int score;              // Coding, Non-coding or Boundary
    private int codingFrame;        // Coding frame
    private int startCodon;         // Position of the first base of codon (WRT Genome)
    private int stopCodon;          // Position of the last base of codon (WRT Genome)
    private String nucSequence;

    public ExonInfo(String transcriptId, String exonId, String chromosome,
        int start, int stop, int direction) {

        this.transcriptId = transcriptId;
        this.exonId       = exonId;
        this.chromosome   = chromosome;
        this.start        = start;
        this.stop         = stop;
        this.direction    = direction;

        this.score        = -100; // Default to arbitrary number...Will use -100
        this.codingFrame  = -100; // Default to arbitrary number...Will use -100
        this.startCodon   = -100;
        this.stopCodon    = -100;
    }

    public ExonInfo(ExonInfo e) {
        this(e.getTranscriptId(), e.getId(), e.getChromosome(),
            e.getStart(), e.getStop(), e.getDirection());
    }

    // Get Functions
    public String getTranscriptId() {
        return this.transcriptId;
    }

    public String getId() {
        return this.exonId;
    }

    public String getChromosome() {
        return this.chromosome;
    }

    public int getStart() {
        return start;
    }

    public int getStop() {
        return stop;
    }

    public int getDirection() {
        return this.direction;
    }

    public String getDirectionStr() {
        if (getDirection() == 1) {
            return GenomeConstant.FORWARD_FLAG;
        }
        return GenomeConstant.REVERSE_FLAG;
    }

    public int getScore() {
        return this.score;
    }

    public int getCodingFrame() {
        return this.codingFrame;
    }

    public int getStartCodon() {
        return this.startCodon;
    }

    public int getStopCodon() {
        return this.stopCodon;
    }

    // Has Functions
    public boolean hasScore() {
        if (score != -100) {
            return true;
        }
        return false;
    }

    public boolean hasCodingFrame() {
        if (codingFrame != -100) {
            return true;
        }
        return false;
    }

    public boolean containsStartCodon() {
        if (startCodon != -100) {
            return true;
        }
        return false;
    }

    public boolean containsStopCodon() {
        if (stopCodon != -100) {
            return true;
        }
        return false;
    }

    public boolean isForward() {
        return getDirection() == 1;
    }

    // Set Functions
    public void setId(String exonId) {
        this.exonId = exonId;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public void setStop(int stop) {
        this.stop = stop;
    }

    public void setScore(int score) {
        this.score = score;
    }

    public void setCodingFrame(int codingFrame) {
        this.codingFrame = codingFrame;
    }

    public void setStartCodon(int startCodon) {
        this.startCodon = startCodon;
    }

    public void setStartCodon(CodonInfo startCodon) {
        setStartCodon(startCodon.getStart());
        if (!isForward()) {
            setStartCodon(startCodon.getStop());
        }
    }

    public void setStopCodon(int stopCodon) {
        this.stopCodon = stopCodon;
    }

    public void setStopCodon(CodonInfo stopCodon) {
        setStopCodon(stopCodon.getStop());
        if (!isForward()) {
            setStopCodon(stopCodon.getStart());
        }
    }

    @Override
    public int hashCode() {
        return new HashCodeBuilder(17, 37).append(chromosome)
                .append(direction)
                .append(start)
                .append(stop).toHashCode();
    }

    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof ExonInfo)) {
            return false;
        }
        if (this == obj) {
            return true;
        }
        ExonInfo other = (ExonInfo) obj;
        return new EqualsBuilder().append(chromosome, other.chromosome)
                .append(direction, other.direction)
                .append(start, other.start)
                .append(stop, other.stop).isEquals();
    }

    @Override
    public String toString() {
        StringBuffer out = new StringBuffer();
        out.append(transcriptId + "\t" + exonId + "\t" + chromosome);
        out.append("\t");
        out.append(getDirectionStr() + "\t" + start + "\t" + stop);

        if (hasScore()) {
            out.append("\t");
            out.append("Score: " + score);
        }

        if (hasCodingFrame()) {
            out.append("\t");
            out.append("Coding Frame: " + codingFrame);
        }

        if (containsStartCodon()) {
            out.append("\t");
            out.append("Start codon: " + startCodon);
        }

        if (containsStopCodon()) {
            out.append("\t");
            out.append("Stop codon: " + stopCodon);
        }

        return out.toString();
    }

}
