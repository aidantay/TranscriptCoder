package au.org.intersect.samifier.inferencer;

import java.io.File;
import java.util.List;
import java.util.Map;

import au.org.intersect.samifier.parser.FastaParser;
import au.org.intersect.samifier.parser.FastaParserImpl;
import au.org.intersect.samifier.parser.TranscriptomeParser;
import au.org.intersect.samifier.parser.TranscriptomeParserImpl;
import au.org.intersect.samifier.domain.CodonTranslationTable;
import au.org.intersect.samifier.domain.MegaExonInfo;
import au.org.intersect.samifier.domain.TranscriptInfo;
import au.org.intersect.samifier.domain.Transcriptome;
import au.org.intersect.samifier.inferencer.ExonFinder;
import au.org.intersect.samifier.inferencer.ExonFinderImpl;
import au.org.intersect.samifier.inferencer.TranscriptGenerator;
import au.org.intersect.samifier.inferencer.TranscriptGeneratorException;
import au.org.intersect.samifier.inferencer.TranscriptGeneratorImpl;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import org.junit.Before;
import org.junit.Test;

public final class TranscriptGeneratorImplTest {

    private StringBuilder genomeString;
    private CodonTranslationTable translationTable;
    private Transcriptome transcriptome;
    private List<TranscriptInfo> refTranscripts;

    @Before
    public void before() throws Exception {
        File t = new File("test/resources/test_transcriptome.gtf");
        TranscriptomeParser parser = new TranscriptomeParserImpl();
        transcriptome  = parser.parseTranscriptomeFile(t);
        refTranscripts = transcriptome.getAllTranscripts();

        File g = new File("test/resources/chr10.fa");
        File c = new File("test/resources/standard_code_translation_table.txt");
        FastaParser fastaParser = new FastaParserImpl(g);
        genomeString     = new StringBuilder(fastaParser.readCode("chr10"));
        translationTable = CodonTranslationTable.parseTableFile(c);
    }

    @Test
    public void testForwardMultiExonTranscript() 
            throws TranscriptGeneratorException {
        TranscriptInfo transcript = transcriptome.getTranscript("ENST00000651811");

        TranscriptGenerator transcriptGenerator   = new TranscriptGeneratorImpl(genomeString, translationTable);
        ExonFinder exonFinder                     = new ExonFinderImpl(refTranscripts);
        MegaExonInfo closestKnownExon = exonFinder.getClosestKnownExon(transcript);
        if (closestKnownExon != null) {
            List<TranscriptInfo> transcripts = transcriptGenerator.inferTranscript(transcript, closestKnownExon);
            assertEquals(transcripts.size(), 1);
            assertEquals(transcripts.get(0).getAllExons().size(), 8);
            assertEquals(transcripts.get(0).getAllCodons().size(), 2);
            assertEquals(transcripts.get(0).getStartCodon().getStart(), 110228224);
            assertEquals(transcripts.get(0).getStopCodon().getStop(), 110284987);
        }
    }

    @Test
    public void testReverseMultiExonTranscript() 
            throws TranscriptGeneratorException {
        TranscriptInfo transcript = transcriptome.getTranscript("ENST00000564130");

        TranscriptGenerator transcriptGenerator   = new TranscriptGeneratorImpl(genomeString, translationTable);
        ExonFinder exonFinder                     = new ExonFinderImpl(refTranscripts);
        MegaExonInfo closestKnownExon = exonFinder.getClosestKnownExon(transcript);
        if (closestKnownExon != null) {
            List<TranscriptInfo> transcripts = transcriptGenerator.inferTranscript(transcript, closestKnownExon);
            assertEquals(transcripts.size(), 1);
            assertEquals(transcripts.get(0).getAllExons().size(), 6);
            assertEquals(transcripts.get(0).getAllCodons().size(), 2);
            assertEquals(transcripts.get(0).getStartCodon().getStop(), 48893);
            assertEquals(transcripts.get(0).getStopCodon().getStart(), 47057);
        }
    }

    @Test
    public void testForwardSingleExonTranscript() 
            throws TranscriptGeneratorException {
        TranscriptInfo transcript = transcriptome.getTranscript("ENST00000441178");

        TranscriptGenerator transcriptGenerator   = new TranscriptGeneratorImpl(genomeString, translationTable);
        ExonFinder exonFinder                     = new ExonFinderImpl(refTranscripts);
        MegaExonInfo closestKnownExon = exonFinder.getClosestKnownExon(transcript);
        if (closestKnownExon != null) {
            List<TranscriptInfo> transcripts = transcriptGenerator.inferTranscript(transcript, closestKnownExon);
            assertEquals(transcripts.size(), 1);
            assertEquals(transcripts.get(0).getAllExons().size(), 3);
            assertEquals(transcripts.get(0).getAllCodons().size(), 2);
            assertEquals(transcripts.get(0).getStartCodon().getStart(), 103245997);
            assertEquals(transcripts.get(0).getStopCodon().getStop(), 103246683);
        }
    }

    @Test
    public void testReverseSingleExonTranscript() 
            throws TranscriptGeneratorException {
        TranscriptInfo transcript = transcriptome.getTranscript("ENST00000371019");

        TranscriptGenerator transcriptGenerator   = new TranscriptGeneratorImpl(genomeString, translationTable);
        ExonFinder exonFinder                     = new ExonFinderImpl(refTranscripts);
        MegaExonInfo closestKnownExon = exonFinder.getClosestKnownExon(transcript);
        if (closestKnownExon != null) {
            List<TranscriptInfo> transcripts = transcriptGenerator.inferTranscript(transcript, closestKnownExon);
            assertEquals(transcripts.size(), 1);
            assertEquals(transcripts.get(0).getAllExons().size(), 3);
            assertEquals(transcripts.get(0).getAllCodons().size(), 2);
            assertEquals(transcripts.get(0).getStartCodon().getStop(), 97334572);
            assertEquals(transcripts.get(0).getStopCodon().getStart(), 97333871);
        }
    }

}
