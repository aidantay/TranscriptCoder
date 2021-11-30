package au.org.intersect.samifier.domain;

import java.util.Comparator;

public class CodonInfoComparator implements Comparator<CodonInfo> {

    @Override
    // -1 = o2 RIGHT OF o1
    //  0 = o2 EQUALS o1
    //  1 = o2 LEFT OF o1   
    public int compare(CodonInfo o1, CodonInfo o2) {
        float o1Id = Float.parseFloat(o1.getExonId());
        float o2Id = Float.parseFloat(o2.getExonId());
        if (o1Id < o2Id) {
            return -1;
        } else if (o1Id > o2Id) {
            return 1;
        } else {
            if (o1.getType().equals(TranscriptomeConstant.START_CODON)
                && o2.getType().equals(TranscriptomeConstant.STOP_CODON)) {
                return -1;
            } else if (o1.getType().equals(TranscriptomeConstant.STOP_CODON)
                && o2.getType().equals(TranscriptomeConstant.START_CODON)) {

            } else {
                if (o1.getType().equals(o2.getType())) {
                    if (o1.getStart() < o2.getStart() && o1.getStop() < o2.getStop()) {
                        return -1;
                    } else if (o1.getStart() > o2.getStart() && o1.getStop() > o2.getStop()) {
                        return 1;
                    }
                }
            }
        }
        return 0;
    }
}
