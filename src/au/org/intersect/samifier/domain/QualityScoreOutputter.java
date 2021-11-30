package au.org.intersect.samifier.domain;

public class QualityScoreOutputter {

    private Transcriptome targetTranscriptome;
    private Transcriptome outputTranscriptome;

    public QualityScoreOutputter(Transcriptome outputTranscriptome, Transcriptome targetTranscriptome) {
            this.outputTranscriptome = outputTranscriptome;
            this.targetTranscriptome = targetTranscriptome;
    }

    public String writeQualityScore () {
    StringBuffer out = new StringBuffer();
    // out.append("TranscriptID\tScore");
    // out.append(System.getProperty("line.separator"));        
    // for (TranscriptInfo transcript : targetTranscriptome.getTranscripts()) {
    //      out.append(transcript.getId());
    //      out.append("\t");
    //      if (outputTranscriptome.hasTranscript(transcript.getId())) {
    //              out.append(outputTranscriptome.getTranscript(transcript.getId()).getQuality());
    //      } else {
    //              out.append("-5");
    //      }
    //     out.append(System.getProperty("line.separator"));
    // }
            return out.toString();
    }

}
