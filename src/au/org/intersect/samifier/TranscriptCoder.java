package au.org.intersect.samifier;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.log4j.FileAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.PatternLayout;

import au.org.intersect.samifier.runner.TranscriptCoderRunner;

public class TranscriptCoder {

    public static void main(String[] args) {      
        // Input Parameters
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Subject/Reference transcript file in gtf format.");
        OptionBuilder.withArgName("Subject transcript File");       
        OptionBuilder.isRequired();
        Option referenceFileOpt = OptionBuilder.create("s");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Query/Target transcript file in gtf format.");
        OptionBuilder.withArgName("Query transcript File");       
        OptionBuilder.isRequired();
        Option targetFileOpt = OptionBuilder.create("q");
        OptionBuilder.hasArg();
        OptionBuilder.withArgName("Translation Table File");
        OptionBuilder.withDescription("File containing a mapping of codons to amino acids, in the format used by NCBI.");
        OptionBuilder.isRequired();
        Option translationTableOpt = OptionBuilder.create("t");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Directory containing the chromosome files in FASTA format for the given genome.");
        OptionBuilder.withArgName("chromosomeDir");
        OptionBuilder.isRequired();
        Option chrDirOpt = OptionBuilder.create("c");

        // Output Parameters
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Database name.");
        OptionBuilder.withArgName("Database Name");
        OptionBuilder.isRequired();
        Option databaseNameOpt = OptionBuilder.create("d");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Filename to write the FASTA format file to.");
        OptionBuilder.withArgName("Output File");
        OptionBuilder.isRequired();
        Option databaseFileOpt = OptionBuilder.create("o");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Filename to write the log into.");
        OptionBuilder.withArgName("Log File");
        OptionBuilder.isRequired();
        Option logFileOpt = OptionBuilder.create("l");

        // Optional Output Parameters
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Filename to write the GFF file to.");
        OptionBuilder.withArgName("GFF File");
        OptionBuilder.isRequired(false);
        Option gffFileOpt = OptionBuilder.create("p");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Filename to write the accession file to.");
        OptionBuilder.withArgName("Accession File");
        OptionBuilder.isRequired(false);
        Option accessionFileOpt = OptionBuilder.create("m");

        Options options = new Options();
        options.addOption(translationTableOpt);
        options.addOption(chrDirOpt);
        options.addOption(databaseNameOpt);
        options.addOption(databaseFileOpt);
        options.addOption(gffFileOpt);
        options.addOption(accessionFileOpt);
        options.addOption(referenceFileOpt);
        options.addOption(targetFileOpt);
        options.addOption(logFileOpt);

        CommandLineParser parser = new GnuParser();
        try {
            CommandLine line = parser.parse(options, args);
            File referenceTranscriptFile = new File(line.getOptionValue("s"));
            File targetTranscriptFile = new File(line.getOptionValue("q"));
            File chromosomeDir = new File(line.getOptionValue("c"));
            File translationTableFile = new File(line.getOptionValue("t"));
            
            String databaseName = line.getOptionValue("d");
            File databaseFile = new File(line.getOptionValue("o"));
            Writer databaseWriter = new FileWriter(databaseFile);

            String logFileName = line.getOptionValue("l");
            setFileLogger(logFileName);

            String gffFilename = line.getOptionValue("p");
            Writer gffWriter = setWriterFile(gffFilename);

            String accessionFilename = line.getOptionValue("m");
            Writer accessionWriter = setWriterFile(accessionFilename);

            TranscriptCoderRunner runner = new TranscriptCoderRunner(
                    referenceTranscriptFile, targetTranscriptFile, 
                    translationTableFile, chromosomeDir,
                    databaseName, databaseWriter, 
                    gffWriter, accessionWriter);

            runner.run();


        } catch (ParseException pe) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("transcriptcoder", options, true);
            System.exit(1);

        } catch (Exception e) {
            System.err.println(e);
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("transcriptcoder", options, true);
            e.printStackTrace();
            System.exit(1);
        }
    }

    private static Writer setWriterFile(String fileName) throws IOException {
        Writer writer = null;
        if (fileName != null) {
            File file = new File(fileName);
            writer = new FileWriter(file);
        }
        
        return writer;
    }

    private static void setFileLogger(String logFileName) {
        Logger.getRootLogger().removeAppender("stdout");
        FileAppender fa = new FileAppender();
        fa.setName("FileLogger");
        fa.setFile(logFileName);
        fa.setLayout(new PatternLayout("%d %-5p %c - %m%n"));
        fa.setThreshold(Level.DEBUG);
        fa.setAppend(true);
        fa.activateOptions();
        Logger.getRootLogger().addAppender(fa);
    }

}
