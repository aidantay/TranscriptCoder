package au.org.intersect.samifier.parser;

import au.org.intersect.samifier.domain.Transcriptome;

import java.io.File;

public interface TranscriptomeParser {

    int CHROMOSOME_PART = 0;
    int SOURCE_PART = 1;
    int TYPE_PART = 2;
    int ATTRIBUTES_PART = 8;
    int START_PART = 3;
    int STOP_PART = 4;
    int STRAND_PART = 6;

    Transcriptome parseTranscriptomeFile(File transcriptomeFile) throws TranscriptomeFileParsingException;
}
