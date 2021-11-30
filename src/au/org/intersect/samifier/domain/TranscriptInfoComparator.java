package au.org.intersect.samifier.domain;

import java.util.Comparator;

public class TranscriptInfoComparator implements Comparator<TranscriptInfo> {

	@Override
	public int compare(TranscriptInfo o1, TranscriptInfo o2) {
		return String.CASE_INSENSITIVE_ORDER.compare(o1.getChromosome(), o2.getChromosome());
	}

	
	
}
