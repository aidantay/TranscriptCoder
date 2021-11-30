package au.org.intersect.samifier.domain;

public class GffOutputterGenerator implements
        ProteinLocationBasedOutputterGenerator {
    private String genomeFilename;

    public GffOutputterGenerator(String genomeFilename) {
        this.genomeFilename = genomeFilename;
    }
    
    public GffOutputterGenerator() {
    	
    }

    @Override
    public GffOutputter getOutputterFor(ProteinLocation proteinLocation) {
        return new GffOutputter(proteinLocation, genomeFilename);
    }

	@Override
	public GffOutputter getOutputterFor(TranscriptInfo transcriptinfo) {
		return new GffOutputter(transcriptinfo);
	}


}
