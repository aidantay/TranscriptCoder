package au.org.intersect.samifier.parser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.log4j.Logger;

import au.org.intersect.samifier.domain.ExonInfo;
import au.org.intersect.samifier.domain.CodonInfo;
import au.org.intersect.samifier.domain.TranscriptInfo;
import au.org.intersect.samifier.domain.Transcriptome;

public class TranscriptomeParserImpl implements TranscriptomeParser {

    private static Logger LOG = Logger.getLogger(TranscriptomeParserImpl.class);

    private static final String CHROMOSOME = "chr";
    private static final Pattern CHROMOSOME_RE = Pattern.compile("^" + CHROMOSOME);

    private static final String STRAND_FORWARD = "+";
    private static final Pattern STRAND_RE = Pattern.compile("^([" + STRAND_FORWARD + "]|[-])$");

    private static final String TRANSCRIPT_TYPE = "transcript";
    private static final Pattern TRANSCRIPT_TYPE_RE = Pattern.compile("^" + TRANSCRIPT_TYPE + "$");

    private static final String EXON_TYPE = "exon";
    private static final Pattern EXON_TYPE_RE = Pattern.compile("^" + EXON_TYPE + "$");

    private static final String CODON_TYPE = "(start_codon|stop_codon)";
    private static final Pattern CODON_TYPE_RE = Pattern.compile("^" + CODON_TYPE + "$");

    private static final String TRANSCRIPT_ID_TAG = "transcript_id";
    private static final Pattern TRANSCRIPT_ID_ATTRIBUTE_RE = Pattern.compile(TRANSCRIPT_ID_TAG + " \"([a-zA-Z0-9\\.\\-_]+)\";");

    private static final String EXON_ID_TAG = "exon_number";
    private static final Pattern EXON_ID_ATTRIBUTE_RE = Pattern.compile(EXON_ID_TAG + " \"([0-9\\-_]+)\";");

    private String transcriptomeFileName;
    private int lineNumber = 0;
    private String line;

    public TranscriptomeParserImpl() {

    }

    public Transcriptome parseTranscriptomeFile(File transcriptomeFile)
            throws TranscriptomeFileParsingException {
        try {
            transcriptomeFileName = transcriptomeFile.getAbsolutePath();
            return doParsing(transcriptomeFile);
        } catch (IOException e) {
            throw new TranscriptomeFileParsingException(e.getMessage());
        }
    }

    // Creates a transcript from the exons and assigns the start and stop codons to the transript
    // Once transcript is 'ready', all exons are merged into the transcriptome (Map of all Exons)
    private Transcriptome doParsing(File transcriptFile)
            throws IOException, TranscriptomeFileParsingException {
      
        Transcriptome transcriptome = new Transcriptome();
        BufferedReader reader = null;
        
        try {
            reader = new BufferedReader(new FileReader(transcriptFile));

            while ((line = reader.readLine()) != null) {
                lineNumber++;
                if (line.startsWith("#")) {
                    continue;
                }

                // Chromosome, source, type, start, stop, score, strand, phase, attributes
                String[] parts = line.split("\\t", 9);
                if (parts.length < 9) {
                    String line1 = ">Line " + lineNumber + ":" + "Not in expected format.\n";
                    String line2 = ">" + line + "\n";
                    String message = line1 + line2;
                    throwParsingException(message);
                }

                // Any exons with no direction specified will throw an exception
                // Users should remove these cases before running the tool
                String type = parts[TYPE_PART];
                if (parseStrand(parts[STRAND_PART]) == 0 || type == null) {
                    String line1 = ">Line " + lineNumber + ":" + "Direction or type not specified.\n";
                    String line2 = ">" + line + "\n";
                    String message = line1 + line2;
                    throwParsingException(message);
                }

                if (TRANSCRIPT_TYPE_RE.matcher(type).matches()) {
                    TranscriptInfo transcript = parseTranscript(parts);
                    processTranscript(transcriptome, transcript);

                } else if (EXON_TYPE_RE.matcher(type).matches()) {
                    ExonInfo exon = parseExon(parts);
                    processExon(transcriptome, exon);

                } else if (CODON_TYPE_RE.matcher(type).matches()) {
                    CodonInfo codon = parseCodon(parts);
                    processCodon(transcriptome, codon);

                } else {
                    LOG.warn(">Line " + lineNumber + ": Skipping unsupported type");
                    continue;
                }
            }

        } finally {
            if (reader != null) {
                reader.close();
            }
        }
        transcriptome.verify();
        return transcriptome;
    }
    
    private void throwParsingException(String message)
            throws TranscriptomeFileParsingException {

        message = "Error in " + transcriptomeFileName + ":\n\t" + message;
        throw new TranscriptomeFileParsingException(message);
    }

    private int parseStrand(String direction) {
        if (!STRAND_RE.matcher(direction).matches()) {
          return 0;
        }
        return STRAND_FORWARD.equals(direction) ? 1 : -1;
    }

    protected TranscriptInfo parseTranscript(String[] parts)
            throws TranscriptomeFileParsingException {

        String chromosome   = extractChromosome(parts[CHROMOSOME_PART]);
        String origin       = parts[SOURCE_PART];
        int start           = Integer.parseInt(parts[START_PART]);
        int stop            = Integer.parseInt(parts[STOP_PART]);
        int direction       = parseStrand(parts[STRAND_PART]);
        String transcriptId = extractAttribute(parts[ATTRIBUTES_PART], TRANSCRIPT_ID_ATTRIBUTE_RE);
        return new TranscriptInfo(transcriptId, chromosome, origin, start, stop, direction);
    }

    protected ExonInfo parseExon(String[] parts)
            throws TranscriptomeFileParsingException {

        String chromosome   = extractChromosome(parts[CHROMOSOME_PART]);
        int start           = Integer.parseInt(parts[START_PART]);
        int stop            = Integer.parseInt(parts[STOP_PART]);
        int direction       = parseStrand(parts[STRAND_PART]);
        String transcriptId = extractAttribute(parts[ATTRIBUTES_PART], TRANSCRIPT_ID_ATTRIBUTE_RE);
        String exonId       = extractAttribute(parts[ATTRIBUTES_PART], EXON_ID_ATTRIBUTE_RE);
        return new ExonInfo(transcriptId, exonId, chromosome, start, stop, direction);
    }

    protected CodonInfo parseCodon(String[] parts)
            throws TranscriptomeFileParsingException {

        String chromosome   = extractChromosome(parts[CHROMOSOME_PART]);
        String type         = parts[TYPE_PART];
        int start           = Integer.parseInt(parts[START_PART]);
        int stop            = Integer.parseInt(parts[STOP_PART]);
        int direction       = parseStrand(parts[STRAND_PART]);
        String transcriptId = extractAttribute(parts[ATTRIBUTES_PART], TRANSCRIPT_ID_ATTRIBUTE_RE);
        String exonId       = extractAttribute(parts[ATTRIBUTES_PART], EXON_ID_ATTRIBUTE_RE);
        return new CodonInfo(transcriptId, exonId, type, chromosome, start, stop, direction);
    }
    
    private void processTranscript(Transcriptome transcriptome, TranscriptInfo transcript) {
        transcriptome.addTranscript(transcript);
    }

    private void processExon(Transcriptome transcriptome, ExonInfo exon)
            throws TranscriptomeFileParsingException {

        TranscriptInfo transcript = transcriptome.getTranscript(exon.getTranscriptId());
        if (transcript.hasExon(exon)) {
            String message = "Transcript " + transcript.getId() 
                + " already contains exon. Unexpected case TODO.";
            throwParsingException(message);
        }
        transcriptome.getTranscript(exon.getTranscriptId()).addExon(exon);
    }

    private void processCodon(Transcriptome transcriptome, CodonInfo codon)
            throws TranscriptomeFileParsingException {

        TranscriptInfo transcript = transcriptome.getTranscript(codon.getTranscriptId());
        if (transcript.hasCodon(codon)) {
            String message = "Transcript " + transcript.getId() 
                + " already contains exon. Unexpected case TODO.";
            throwParsingException(message);
        }
        transcriptome.getTranscript(codon.getTranscriptId()).addCodon(codon);
    }
    
    private String extractChromosome(String attributes) {
        Matcher m = CHROMOSOME_RE.matcher(attributes);
        if (m.find()) {
            return attributes;
        }
        return "chr" + attributes; // make compiler happy
    }

    private String extractAttribute(String attributes, Pattern p)
            throws TranscriptomeFileParsingException {

        Matcher m = p.matcher(attributes);
        if (m.find()) {
            return m.group(1);
        }
        throwParsingException("Attribute transcript_id not found");
        return null; // make compiler happy
    }
}
