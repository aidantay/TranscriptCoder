package au.org.intersect.samifier.domain;

public class AccessionOutputter implements Outputter {
    private String name;

    public AccessionOutputter(ProteinLocation proteinLocation) {
        this.name = proteinLocation.getName();
    }
    
    public String getOutput() throws OutputException {
        return name + " " + name + " " + name
                + System.getProperty("line.separator");
    }
    
    // For generating the protein sequence of a transcript isoform
    public AccessionOutputter(TranscriptInfo transcriptInfo) {
    	this.name = transcriptInfo.getId();
    }

	public String getTranscriptOutput() throws OutputException {
		return name + " " + name + " " + name
				+ System.getProperty("line.separator");
    }
    
}
