package au.org.intersect.samifier.domain;

public class CodonInfo {

    private String transcriptId;
    private String exonId;
    private String type;
    private String chromosome;
    private int start;              // Start position of the transcript (WRT Genome)
    private int stop;               // End position of the transcript   (WRT Genome)
    private int direction;

    public CodonInfo(String transcriptId, String exonId, String type,
        String chromosome, int start, int stop, int direction) {

        this.transcriptId = transcriptId;
        this.exonId       = exonId;
        this.type         = type;
        this.chromosome   = chromosome;
        this.start        = start;
        this.stop         = stop;
        this.direction    = direction;
    }

    // Get Functions   
    public String getTranscriptId() {
        return this.transcriptId;
    }

    public String getExonId() {
        return this.exonId;
    }

    public String getType() {
        return this.type;
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
        return direction;
    }

    public String getDirectionStr() {
        if (getDirection() == 1) {
            return GenomeConstant.FORWARD_FLAG;
        }
        return GenomeConstant.REVERSE_FLAG;
    }

    public boolean isStart() {
        if (getType().equals(TranscriptomeConstant.START_CODON)) {
            return true;
        }
        return false;
    }

    public boolean isStop() {
        if (getType().equals(TranscriptomeConstant.STOP_CODON)) {
            return true;
        }
        return false;
    }

    @Override
    public String toString() {
        StringBuffer out = new StringBuffer();
        out.append(transcriptId + "\t" + exonId + "\t" + type + "\t" + chromosome);
        out.append("\t");
        out.append(getDirectionStr() + "\t" + start + "\t" + stop);
        return out.toString();
    }
}
