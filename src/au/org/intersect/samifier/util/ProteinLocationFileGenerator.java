package au.org.intersect.samifier.util;

import au.org.intersect.samifier.domain.OutputException;
import au.org.intersect.samifier.domain.ProteinLocation;
import au.org.intersect.samifier.domain.ProteinLocationBasedOutputterGenerator;
import au.org.intersect.samifier.domain.TranscriptInfo;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.List;

public class ProteinLocationFileGenerator {
    public static void generateFile(List<ProteinLocation> locations,
        Writer fileWriter, ProteinLocationBasedOutputterGenerator outputterGenerator)
            throws IOException {
        generateFile(locations, fileWriter, outputterGenerator, null);
    }

    public static void generateFile(List<ProteinLocation> locations,
        Writer fileWriter, ProteinLocationBasedOutputterGenerator outputterGenerator,
        String header)
            throws IOException {

        if (fileWriter != null) {
            BufferedWriter writer = null;
            try {
                writer = new BufferedWriter(fileWriter);
                if (header != null && header.length() > 0) {
                    writer.append(header);
                    writer.append("\n");
                }
                for (ProteinLocation location : locations) {
                    writer.append(outputterGenerator.getOutputterFor(location)
                            .getOutput());
                }
            } catch (IOException e) {
                System.exit(-1);
            } catch (OutputException e) {
                System.err.println(e);
                System.exit(-1);
            } finally {
                if (writer != null) {
                    writer.flush();
                }
            }
        }
    }
    
    // For generating the protein sequence of a transcript isoform
    public static void generateTranscriptSequenceFile(List<TranscriptInfo> transcripts,
        Writer fileWriter, ProteinLocationBasedOutputterGenerator outputterGenerator)
            throws IOException, OutputException {
        generateTranscriptSequenceFile(transcripts, fileWriter, outputterGenerator, null);
    }
    
    public static void generateTranscriptSequenceFile(List<TranscriptInfo> transcripts,
        Writer fileWriter, ProteinLocationBasedOutputterGenerator outputterGenerator,
        String header)
            throws IOException, OutputException {

        if (fileWriter != null) {
            BufferedWriter writer = null;
            try {
                writer = new BufferedWriter(fileWriter);
                if (header != null && header.length() > 0) {
                    writer.append(header);
                    writer.append("\n");
                }
                for (TranscriptInfo transcriptinfo : transcripts) {
                    writer.append(outputterGenerator.getOutputterFor(transcriptinfo).getTranscriptOutput());
                }
            } catch (IOException e) {
                System.exit(-1);
            } catch (OutputException e) {
                System.err.println(e);
                System.exit(-1);
            } finally {
                if (writer != null) {
                    writer.flush();
                }
            }
        }
    }
}
