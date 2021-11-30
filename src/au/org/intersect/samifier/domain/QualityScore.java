package au.org.intersect.samifier.domain;

public class QualityScore {

	private int exonScore;
	private int phaseScore;
	private int codonPairScore;
	
	// Set functions
	public QualityScore() {
		this.exonScore = 0;
		this.phaseScore = 0;
		this.codonPairScore = 0;
	}

	public void setExonScore(int exonScore) {
		this.exonScore = exonScore;
	}

	public void setPhaseScore(int phaseScore) {
		this.phaseScore = phaseScore;
	}

	public void setCodonPairScore(int codonPairScore) {
		this.codonPairScore = codonPairScore;
	}

	// Get functions
	public int getExonScore() {
		return exonScore;
	}

	public int getPhaseScore() {
		return phaseScore;
	}

	public int getCodonPairScore() {
		return codonPairScore;
	}

	public int getTotal() {
		return exonScore - phaseScore - codonPairScore;
	}

	public String toString() {
		StringBuffer out = new StringBuffer();
		out.append(exonScore + "\t" + phaseScore + "\t" + codonPairScore + "\t" + getTotal());
		return out.toString();
	}
	
}
