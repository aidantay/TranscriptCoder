package au.org.intersect.samifier.domain;

import java.util.Comparator;

public class ExonInfoComparator implements Comparator<ExonInfo> {

    @Override
    // -1 = o2 RIGHT OF o1
    //  0 = o2 EQUALS o1
    //  1 = o2 LEFT OF o1   
    public int compare(ExonInfo o1, ExonInfo o2) {
        float o1Id = Float.parseFloat(o1.getId());
        float o2Id = Float.parseFloat(o2.getId());
        if (o1Id < o2Id) {
            return -1;
        } else if (o1Id > o2Id) {
            return 1;
        } else {
            if (o1.getStart() < o2.getStart() && o1.getStop() < o2.getStop()) {
                return -1;
            } else if (o1.getStart() > o2.getStart() && o1.getStop() > o2.getStop()) {
                return 1;
            }            
        }        
        return 0;
    }
}
