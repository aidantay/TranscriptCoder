package au.org.intersect.samifier.domain;

/**
 * Created with IntelliJ IDEA. User: diego Date: 25/10/12 Time: 10:06 AM To
 * change this template use File | Settings | File Templates.
 */
public class AccessionOutputterGenerator implements
        ProteinLocationBasedOutputterGenerator {
    public AccessionOutputterGenerator() {
    }

    public Outputter getOutputterFor(ProteinLocation proteinLocation) {
        return new AccessionOutputter(proteinLocation);
    }

	public Outputter getOutputterFor(TranscriptInfo transcriptInfo) {
		return new AccessionOutputter(transcriptInfo);
	}

}
