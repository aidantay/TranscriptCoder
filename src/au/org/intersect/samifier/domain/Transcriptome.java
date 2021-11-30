package au.org.intersect.samifier.domain;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import org.apache.commons.lang3.tuple.Pair; 

import au.org.intersect.samifier.parser.TranscriptomeFileParsingException;

public class Transcriptome {

    private static Logger LOG = Logger.getLogger(Transcriptome.class);

    // Transcript Table
    // TranscriptId|Chromosome|Strand|Start|End
    // ENST1       |1         |+     |1    |4
    // ENST2       |2         |-     |5    |10
    private Map<String, TranscriptInfo> transcriptMap;

    public Transcriptome() {
        this.transcriptMap = new HashMap<String, TranscriptInfo>();
    }

    // Get Functions        
    public TranscriptInfo getTranscript(String transcriptId) {
        return transcriptMap.get(transcriptId);
    }

    public List<TranscriptInfo> getAllTranscripts() {
        List<TranscriptInfo> transcripts = new ArrayList<TranscriptInfo>(transcriptMap.values());
        Collections.sort(transcripts, new TranscriptInfoComparator());
        return transcripts;
    }

    public List<TranscriptInfo> getAllTranscripts(String chromosome) {
        List<TranscriptInfo> transcripts = new ArrayList<TranscriptInfo>();        
        for (TranscriptInfo t : getAllTranscripts()) {
            if (t.getChromosome().equals(chromosome)) {
                transcripts.add(t);
            }
        }
        return transcripts;
    }

    // Add Functions
    public void addTranscript(TranscriptInfo transcript) {
        transcriptMap.put(transcript.getId(), transcript);
    }

    // Has Functions
    public boolean hasTranscript(String transcriptId) {
        return transcriptMap.containsKey(transcriptId);
    }

    public void verify()
            throws TranscriptomeFileParsingException {

        Iterator<String> iter = transcriptMap.keySet().iterator();
        while (iter.hasNext()) {
            TranscriptInfo t = getTranscript(iter.next());
            // Transcripts that are missing a codon will be removed
            if (!t.isValid()) {
                LOG.warn("Transcript " + t.getId() + " has a missing codon (either start or stop). Removing transcript.");
                iter.remove();

            } else {
                t.update();
            }
        }
    }

    public Map<String, List<TranscriptInfo>> sort() {
        Map<String, List<TranscriptInfo>> locationByChromosome = new HashMap<String, List<TranscriptInfo>>();

        for (TranscriptInfo t : getAllTranscripts()) {
            String chromosome = t.getChromosome();
            List<TranscriptInfo> list = locationByChromosome.get(chromosome);
            if (list == null) {
                list = new ArrayList<TranscriptInfo>();
            }
            list.add(t);
            locationByChromosome.put(chromosome, list);
        }
        return locationByChromosome;
    }

    public String toString() {
        StringBuffer out = new StringBuffer();

        for (String transcriptId : transcriptMap.keySet()) {
            TranscriptInfo transcript = getTranscript(transcriptId);
            out.append(transcript);
            out.append(System.getProperty("line.separator"));

            List<ExonInfo> exons = transcript.getAllExons();
            for (ExonInfo e : exons) {
                out.append("\t");
                out.append(e);
                out.append(System.getProperty("line.separator"));
            }

            List<CodonInfo> codons = transcript.getAllCodons();
            for (CodonInfo c : codons) {
                out.append("\t\t");
                out.append(c);                   
                out.append(System.getProperty("line.separator"));
            }
        }
        return out.toString();
    }

}
