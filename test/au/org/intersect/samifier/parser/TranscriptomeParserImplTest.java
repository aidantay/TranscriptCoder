package au.org.intersect.samifier.parser;

import java.io.File;

import au.org.intersect.samifier.domain.TranscriptInfo;
import au.org.intersect.samifier.domain.Transcriptome;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Before;
import org.junit.Test;

public final class TranscriptomeParserImplTest {

    private File f;

    @Before
    public void before() throws Exception {
        f = new File("test/resources/test_transcriptome.gtf");
    }
    
    @Test
    public void testParseTranscriptomeFile()
            throws TranscriptomeFileParsingException {

        TranscriptomeParser parser = new TranscriptomeParserImpl();
        Transcriptome transcriptome = parser.parseTranscriptomeFile(f);

        assertTrue(transcriptome.hasTranscript("ENST00000651811"));
        assertTrue(transcriptome.hasTranscript("ENST00000441178"));
        assertTrue(transcriptome.hasTranscript("ENST00000416477"));
        assertTrue(transcriptome.hasTranscript("ENST00000564130"));
        assertTrue(transcriptome.hasTranscript("ENST00000371019"));
        assertTrue(transcriptome.hasTranscript("ENST00000562162"));
    }

    @Test
    public void testForwardMultiExonCodingTranscript()
            throws TranscriptomeFileParsingException {

        TranscriptomeParser parser = new TranscriptomeParserImpl();
        Transcriptome transcriptome = parser.parseTranscriptomeFile(f);

        TranscriptInfo t = transcriptome.getTranscript("ENST00000651811");
        assertEquals(t.getStart(), 110226817);
        assertEquals(t.getStop(), 110286585);
        assertEquals(t.getDirectionStr(), "+");
        assertEquals(t.getAllExons().size(), 6);
        assertEquals(t.getAllCodons().size(), 2);
        assertEquals(t.getStartCodon().getStart(), 110228224);
        assertEquals(t.getStopCodon().getStop(), 110284987);
    }

    @Test
    public void testForwardSingleExonCodingTranscript()
            throws TranscriptomeFileParsingException {

        TranscriptomeParser parser = new TranscriptomeParserImpl();
        Transcriptome transcriptome = parser.parseTranscriptomeFile(f);

        TranscriptInfo t = transcriptome.getTranscript("ENST00000441178");
        assertEquals(t.getStart(), 103245887);
        assertEquals(t.getStop(), 103248016);
        assertEquals(t.getDirectionStr(), "+");
        assertEquals(t.getAllExons().size(), 1);
        assertEquals(t.getAllCodons().size(), 2);
        assertEquals(t.getStartCodon().getStart(), 103245997);
        assertEquals(t.getStopCodon().getStop(), 103246683);
    }

    @Test
    public void testForwardMultiExonNoncodingTranscript()
            throws TranscriptomeFileParsingException {

        TranscriptomeParser parser = new TranscriptomeParserImpl();
        Transcriptome transcriptome = parser.parseTranscriptomeFile(f);

        TranscriptInfo t = transcriptome.getTranscript("ENST00000416477");
        assertEquals(t.getStart(), 44712);
        assertEquals(t.getStop(), 46884);
        assertEquals(t.getDirectionStr(), "+");
        assertEquals(t.getAllExons().size(), 6);
        assertEquals(t.getAllCodons().size(), 0);
    }

    @Test
    public void testReverseMultiExonCodingTranscript()
            throws TranscriptomeFileParsingException {

        TranscriptomeParser parser = new TranscriptomeParserImpl();
        Transcriptome transcriptome = parser.parseTranscriptomeFile(f);

        TranscriptInfo t = transcriptome.getTranscript("ENST00000564130");
        assertEquals(t.getStart(), 46892);
        assertEquals(t.getStop(), 74163);
        assertEquals(t.getDirectionStr(), "-");
        assertEquals(t.getAllExons().size(), 4);
        assertEquals(t.getAllCodons().size(), 2);
        assertEquals(t.getStartCodon().getStop(), 48893);
        assertEquals(t.getStopCodon().getStart(), 47057);
    }

    @Test
    public void testReverseSingleExonCodingTranscript()
            throws TranscriptomeFileParsingException {

        TranscriptomeParser parser = new TranscriptomeParserImpl();
        Transcriptome transcriptome = parser.parseTranscriptomeFile(f);

        TranscriptInfo t = transcriptome.getTranscript("ENST00000371019");
        assertEquals(t.getStart(), 97332497);
        assertEquals(t.getStop(), 97334729);
        assertEquals(t.getDirectionStr(), "-");
        assertEquals(t.getAllExons().size(), 1);
        assertEquals(t.getAllCodons().size(), 2);
        assertEquals(t.getStartCodon().getStop(), 97334572);
        assertEquals(t.getStopCodon().getStart(), 97333871);
    }

    @Test
    public void testReverseMultiExonNoncodingTranscript()
            throws TranscriptomeFileParsingException {

        TranscriptomeParser parser = new TranscriptomeParserImpl();
        Transcriptome transcriptome = parser.parseTranscriptomeFile(f);

        TranscriptInfo t = transcriptome.getTranscript("ENST00000562162");
        assertEquals(t.getStart(), 14061);
        assertEquals(t.getStop(), 14604);
        assertEquals(t.getDirectionStr(), "-");
        assertEquals(t.getAllExons().size(), 2);
        assertEquals(t.getAllCodons().size(), 0);
    }

}
