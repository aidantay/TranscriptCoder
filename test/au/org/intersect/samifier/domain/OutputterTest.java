package au.org.intersect.samifier.domain;

import java.io.IOException;
import java.io.File;
import java.util.List;
import java.util.Map;

import au.org.intersect.samifier.parser.FastaParser;
import au.org.intersect.samifier.parser.FastaParserImpl;
import au.org.intersect.samifier.parser.TranscriptomeParser;
import au.org.intersect.samifier.parser.TranscriptomeParserImpl;
import au.org.intersect.samifier.domain.AccessionOutputterGenerator;
import au.org.intersect.samifier.domain.CodonTranslationTable;
import au.org.intersect.samifier.domain.GffOutputterGenerator;
import au.org.intersect.samifier.domain.MegaExonInfo;
import au.org.intersect.samifier.domain.ProteinOutputterGenerator;
import au.org.intersect.samifier.domain.TranscriptInfo;
import au.org.intersect.samifier.domain.Transcriptome;
import au.org.intersect.samifier.inferencer.ExonFinder;
import au.org.intersect.samifier.inferencer.ExonFinderImpl;
import au.org.intersect.samifier.inferencer.TranscriptGenerator;
import au.org.intersect.samifier.inferencer.TranscriptGeneratorImpl;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Before;
import org.junit.Test;

public final class OutputterTest {
    
    private StringBuilder genomeString;
    private CodonTranslationTable translationTable;
    private Transcriptome transcriptome;
    private TranscriptGenerator transcriptGenerator;
    private ExonFinder exonFinder;

    @Before
    public void before() throws Exception {
        File t = new File("test/resources/test_transcriptome.gtf");
        TranscriptomeParser parser = new TranscriptomeParserImpl();
        transcriptome  = parser.parseTranscriptomeFile(t);
        List<TranscriptInfo> refTranscripts = transcriptome.getAllTranscripts();

        File g = new File("test/resources/chr10.fa");
        File c = new File("test/resources/standard_code_translation_table.txt");
        FastaParser fastaParser = new FastaParserImpl(g);
        genomeString     = new StringBuilder(fastaParser.readCode("chr10"));
        translationTable = CodonTranslationTable.parseTableFile(c);

        transcriptGenerator = new TranscriptGeneratorImpl(genomeString, translationTable);
        exonFinder          = new ExonFinderImpl(refTranscripts);
    }

    @Test
    public void testRandom()  
            throws Exception {

        /**********
         * Additional Test cases that we manual inspect
         *
         * ENST00000336378 Forward MultiExon
         * ENST00000463345 Forward MultiExon
         * ENST00000361320 Forward MultiExon
         * ENST00000224950 Reverse MultiExon
         * ENST00000528354 Reverse MultiExon
         * ENST00000361310 Reverse MultiExon
         * ENST00000374115 Forward SingleExon
         * ENST00000520547 Forward SingleExon
         * ENST00000331173 Reverse SingleExon
         * ENST00000371852 Reverse SingleExon
         * ENST00000374314 Reverse SingleExon
         * ENST00000374127 Reverse SingleExon
         * ENST00000399775 Reverse SingleExon
         * ENST00000557985 Reverse SingleExon
         * ENST00000615790 Incomplete
         **********/

        TranscriptInfo transcript = transcriptome.getTranscript("ENST00000615790");
        // TranscriptInfo transcript = null;
        if (transcript != null) {
            MegaExonInfo closestKnownExon = exonFinder.getClosestKnownExon(transcript);
            System.out.println(transcript);
            System.out.println(closestKnownExon);
            System.out.println();
            if (closestKnownExon != null) {
                List<TranscriptInfo> transcripts = transcriptGenerator.inferTranscript(transcript, closestKnownExon);
                System.out.println(transcripts);

                ProteinOutputterGenerator fastaGenerator = new ProteinOutputterGenerator("test", genomeString, translationTable);
                GffOutputterGenerator gffGenerator       = new GffOutputterGenerator();
                AccessionOutputterGenerator idGenerator  = new AccessionOutputterGenerator();
                for (TranscriptInfo t : transcripts) {
                    String fastaLine = fastaGenerator.getOutputterFor(t).getTranscriptOutput();
                    System.out.println(fastaLine);
                    String gffLine = gffGenerator.getOutputterFor(t).getTranscriptOutput();
                    String idLine = idGenerator.getOutputterFor(t).getTranscriptOutput();
                }
            }
        }

    }

    @Test
    public void testForwardMultiExonTranscript()  
            throws Exception {

        TranscriptInfo transcript = transcriptome.getTranscript("ENST00000651811");
        MegaExonInfo closestKnownExon = exonFinder.getClosestKnownExon(transcript);
        if (closestKnownExon != null) {
            List<TranscriptInfo> transcripts = transcriptGenerator.inferTranscript(transcript, closestKnownExon);

            ProteinOutputterGenerator fastaGenerator = new ProteinOutputterGenerator("test", genomeString, translationTable);
            GffOutputterGenerator gffGenerator       = new GffOutputterGenerator();
            AccessionOutputterGenerator idGenerator  = new AccessionOutputterGenerator();
            for (TranscriptInfo t : transcripts) {
                String fastaLine = fastaGenerator.getOutputterFor(t).getTranscriptOutput();
                boolean containsFastaString = fastaLine.contains("MPSPRLQHSKPPRRLSRAQKHSSGSSNTST");
                assertTrue(containsFastaString);

                String gffLine = gffGenerator.getOutputterFor(t).getTranscriptOutput();
                boolean containsGffString = gffLine.contains("TranscriptCoder\tgene\t110226817\t110286585\t0\t+");
                assertTrue(containsGffString);

                String idLine = idGenerator.getOutputterFor(t).getTranscriptOutput();
                boolean containsIdString = idLine.contains("ENST00000651811");
                assertTrue(containsIdString);
            }
        }
    }

    @Test
    public void testReverseMultiExonTranscript()  
            throws Exception {

        TranscriptInfo transcript = transcriptome.getTranscript("ENST00000564130");
        MegaExonInfo closestKnownExon = exonFinder.getClosestKnownExon(transcript);
        if (closestKnownExon != null) {
            List<TranscriptInfo> transcripts = transcriptGenerator.inferTranscript(transcript, closestKnownExon);

            ProteinOutputterGenerator fastaGenerator = new ProteinOutputterGenerator("test", genomeString, translationTable);
            GffOutputterGenerator gffGenerator       = new GffOutputterGenerator();
            AccessionOutputterGenerator idGenerator  = new AccessionOutputterGenerator();
            for (TranscriptInfo t : transcripts) {
                String fastaLine = fastaGenerator.getOutputterFor(t).getTranscriptOutput();
                boolean containsFastaString = fastaLine.contains("MNMPSTPLAPTTGTATCSWSASTCTTTRPA");
                assertTrue(containsFastaString);

                String gffLine = gffGenerator.getOutputterFor(t).getTranscriptOutput();
                boolean containsGffString = gffLine.contains("TranscriptCoder\tgene\t46892\t74163\t0\t-");
                assertTrue(containsGffString);

                String idLine = idGenerator.getOutputterFor(t).getTranscriptOutput();
                boolean containsIdString = idLine.contains("ENST00000564130");
                assertTrue(containsIdString);
            }
        }
    }

    @Test
    public void testForwardSingleExonTranscript()  
            throws Exception {

        TranscriptInfo transcript = transcriptome.getTranscript("ENST00000441178");
        MegaExonInfo closestKnownExon = exonFinder.getClosestKnownExon(transcript);
        if (closestKnownExon != null) {
            List<TranscriptInfo> transcripts = transcriptGenerator.inferTranscript(transcript, closestKnownExon);

            ProteinOutputterGenerator fastaGenerator = new ProteinOutputterGenerator("test", genomeString, translationTable);
            GffOutputterGenerator gffGenerator       = new GffOutputterGenerator();
            AccessionOutputterGenerator idGenerator  = new AccessionOutputterGenerator();
            for (TranscriptInfo t : transcripts) {
                String fastaLine = fastaGenerator.getOutputterFor(t).getTranscriptOutput();
                boolean containsFastaString = fastaLine.contains("MASGCKIGPSILNSDLANLGAKCLQMLDSG");
                assertTrue(containsFastaString);

                String gffLine = gffGenerator.getOutputterFor(t).getTranscriptOutput();
                boolean containsGffString = gffLine.contains("TranscriptCoder\tgene\t103245887\t103248016\t0\t+");
                assertTrue(containsGffString);

                String idLine = idGenerator.getOutputterFor(t).getTranscriptOutput();
                boolean containsIdString = idLine.contains("ENST00000441178");
                assertTrue(containsIdString);
            }
        }
    }

    @Test
    public void testReverseSingleExonTranscript()  
            throws Exception {

        TranscriptInfo transcript = transcriptome.getTranscript("ENST00000371019");
        MegaExonInfo closestKnownExon = exonFinder.getClosestKnownExon(transcript);
        if (closestKnownExon != null) {
            List<TranscriptInfo> transcripts = transcriptGenerator.inferTranscript(transcript, closestKnownExon);

            ProteinOutputterGenerator fastaGenerator = new ProteinOutputterGenerator("test", genomeString, translationTable);
            GffOutputterGenerator gffGenerator       = new GffOutputterGenerator();
            AccessionOutputterGenerator idGenerator  = new AccessionOutputterGenerator();
            for (TranscriptInfo t : transcripts) {
                String fastaLine = fastaGenerator.getOutputterFor(t).getTranscriptOutput();
                boolean containsFastaString = fastaLine.contains("MPCRREEEEEAGEEAEGEEEEDDSFLLLQQ");
                assertTrue(containsFastaString);

                String gffLine = gffGenerator.getOutputterFor(t).getTranscriptOutput();
                boolean containsGffString = gffLine.contains("TranscriptCoder\tgene\t97332497\t97334729\t0\t-");
                assertTrue(containsGffString);

                String idLine = idGenerator.getOutputterFor(t).getTranscriptOutput();
                boolean containsIdString = idLine.contains("ENST00000371019");
                assertTrue(containsIdString);
            }
        }
    }

    @Test
    public void testReverseSingleExonNoUtrTranscript()  
            throws Exception {

        TranscriptInfo transcript = transcriptome.getTranscript("ENST00000374127");
        MegaExonInfo closestKnownExon = exonFinder.getClosestKnownExon(transcript);
        if (closestKnownExon != null) {
            List<TranscriptInfo> transcripts = transcriptGenerator.inferTranscript(transcript, closestKnownExon);

            ProteinOutputterGenerator fastaGenerator = new ProteinOutputterGenerator("test", genomeString, translationTable);
            GffOutputterGenerator gffGenerator       = new GffOutputterGenerator();
            AccessionOutputterGenerator idGenerator  = new AccessionOutputterGenerator();
            for (TranscriptInfo t : transcripts) {
                String fastaLine = fastaGenerator.getOutputterFor(t).getTranscriptOutput();
                boolean containsFastaString = fastaLine.contains("MPRTLSLHEITDLLETDDSIEASAIVIQPP");
                assertTrue(containsFastaString);

                String gffLine = gffGenerator.getOutputterFor(t).getTranscriptOutput();
                boolean containsGffString = gffLine.contains("TranscriptCoder\tgene\t49515333\t49517114\t0\t-");
                assertTrue(containsGffString);

                String idLine = idGenerator.getOutputterFor(t).getTranscriptOutput();
                boolean containsIdString = idLine.contains("ENST00000374127");
                assertTrue(containsIdString);
            }
        }
    }


}
