package au.org.intersect.samifier.inferencer;

import java.util.List;

import au.org.intersect.samifier.domain.MegaExonInfo;
import au.org.intersect.samifier.domain.TranscriptInfo;

public interface TranscriptGenerator {

    public List<TranscriptInfo> inferTranscript(TranscriptInfo transcript, 
        MegaExonInfo closestKnownExon) throws TranscriptGeneratorException;

}
