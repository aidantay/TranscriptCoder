package au.org.intersect.samifier.inferencer;

import java.util.List;

import au.org.intersect.samifier.domain.MegaExonInfo;
import au.org.intersect.samifier.domain.TranscriptInfo;

public interface ExonFinder {

    public MegaExonInfo getClosestKnownExon(TranscriptInfo transcript);

}
