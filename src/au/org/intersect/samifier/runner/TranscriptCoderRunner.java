package au.org.intersect.samifier.runner;

import java.io.File;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import au.org.intersect.samifier.domain.AccessionOutputterGenerator;
import au.org.intersect.samifier.domain.CodonTranslationTable;
import au.org.intersect.samifier.domain.GffOutputterGenerator;
import au.org.intersect.samifier.domain.MegaExonInfo;
import au.org.intersect.samifier.domain.OutputException;
import au.org.intersect.samifier.domain.ProteinOutputterGenerator;
import au.org.intersect.samifier.domain.TranscriptInfo;
import au.org.intersect.samifier.domain.TranscriptInfoComparator;
import au.org.intersect.samifier.domain.Transcriptome;
import au.org.intersect.samifier.domain.TranslationTableParsingException;
import au.org.intersect.samifier.domain.UnknownCodonException;
import au.org.intersect.samifier.parser.FastaParser;
import au.org.intersect.samifier.parser.FastaParserException;
import au.org.intersect.samifier.parser.FastaParserImpl;
import au.org.intersect.samifier.parser.TranscriptomeParserImpl;
import au.org.intersect.samifier.inferencer.ExonFinder;
import au.org.intersect.samifier.inferencer.ExonFinderImpl;
import au.org.intersect.samifier.inferencer.TranscriptGenerator;
import au.org.intersect.samifier.inferencer.TranscriptGeneratorImpl;
import au.org.intersect.samifier.util.ProteinLocationFileGenerator;

public class TranscriptCoderRunner {

    // Command line options
    private File refTranscriptomeFile;
    private File targetTranscriptomeFile;
    private FastaParser fastaParser;
    private CodonTranslationTable translationTable;
    private String databaseName;
    private Writer databaseWriter;
    private Writer gffWriter;
    private Writer accessionWriter;

    private Transcriptome outputTranscriptome;

    public TranscriptCoderRunner (File refTranscriptomeFile, File targetTranscriptomeFile, 
        File translationTableFile, File chromosomeDir, 
        String databaseName, Writer databaseWriter, 
        Writer gffWriter, Writer accessionWriter)
            throws Exception {

        this.refTranscriptomeFile    = refTranscriptomeFile;
        this.targetTranscriptomeFile = targetTranscriptomeFile;
        this.translationTable        = CodonTranslationTable.parseTableFile(translationTableFile);
        this.fastaParser             = new FastaParserImpl(chromosomeDir);
        this.databaseName            = databaseName;
        this.databaseWriter          = databaseWriter;
        this.gffWriter               = gffWriter;
        this.accessionWriter         = accessionWriter;
    }

    public void run()
            throws Exception {

        System.out.print("Parsing...");
        TranscriptomeParserImpl refTranscriptParser = new TranscriptomeParserImpl();
        Transcriptome refTranscriptome              = refTranscriptParser.parseTranscriptomeFile(refTranscriptomeFile);

        TranscriptomeParserImpl targetTranscriptParser = new TranscriptomeParserImpl();
        Transcriptome targetTranscriptome              = targetTranscriptParser.parseTranscriptomeFile(targetTranscriptomeFile);
        System.out.println("Done!");

        System.out.print("Analysing transcriptomes...");
        Map<String, List<TranscriptInfo>> sortedTargetTranscriptome = targetTranscriptome.sort();
        Map<String, List<TranscriptInfo>> sortedRefTranscriptome    = refTranscriptome.sort();
        List<TranscriptInfo> outputTranscripts                      = new ArrayList<TranscriptInfo>();

        for (String chromosome : sortedTargetTranscriptome.keySet()) {
            System.out.println("Transcripts on chromosome " + chromosome);
            List<TranscriptInfo> refTranscripts    = sortedRefTranscriptome.get(chromosome);
            List<TranscriptInfo> targetTranscripts = sortedTargetTranscriptome.get(chromosome);
            List<TranscriptInfo> outputList        = new ArrayList<TranscriptInfo>();

            StringBuilder genomeString = readGenomeFile(chromosome);
            TranscriptGenerator transcriptGenerator = new TranscriptGeneratorImpl(genomeString, translationTable);
            ExonFinder exonFinder                   = new ExonFinderImpl(refTranscripts);

            for (TranscriptInfo prevTranscript : targetTranscripts) {
                MegaExonInfo closestKnownExon = exonFinder.getClosestKnownExon(prevTranscript);
                if (closestKnownExon != null) {
                    List<TranscriptInfo> currTranscripts = transcriptGenerator.inferTranscript(prevTranscript, closestKnownExon);
                    outputList.addAll(currTranscripts);
                }
            }
            generateSequenceDatabase(outputList, genomeString);
            outputTranscripts.addAll(outputList);
        }
        databaseWriter.close();
        System.out.println("Done!");

        generateGffFile(outputTranscripts);
        generateAccessionFile(outputTranscripts);
    }

    private void generateSequenceDatabase(List<TranscriptInfo> outputList, StringBuilder genomeString)
            throws IOException, FastaParserException, TranslationTableParsingException, OutputException {
        
        ProteinOutputterGenerator outputterGenerator = new ProteinOutputterGenerator(databaseName, genomeString, translationTable);
        ProteinLocationFileGenerator.generateTranscriptSequenceFile(outputList, databaseWriter, outputterGenerator);
    }
    
    private void generateGffFile(List<TranscriptInfo> outputTranscripts)
            throws IOException, OutputException {

        if (gffWriter != null) {
            System.out.print("Writing GFF...");
            GffOutputterGenerator outputterGenerator = new GffOutputterGenerator();
            ProteinLocationFileGenerator.generateTranscriptSequenceFile(outputTranscripts, gffWriter,
                    outputterGenerator, "##gff-version 3");
            gffWriter.close();
            System.out.println("Done!");
        }
    }

    private void generateAccessionFile(List<TranscriptInfo> outputTranscripts)
            throws IOException, OutputException {

        if (accessionWriter != null) {
            System.out.print("Wrting Accession file...");
            AccessionOutputterGenerator outputterGenerator = new AccessionOutputterGenerator();
            ProteinLocationFileGenerator.generateTranscriptSequenceFile(outputTranscripts, accessionWriter,
                    outputterGenerator);
            accessionWriter.close();
            System.out.println("Done!");
        }
    }

    private StringBuilder readGenomeFile(String chromosome)
            throws IOException, FastaParserException{

        return new StringBuilder(fastaParser.readCode(chromosome));
    }

}
