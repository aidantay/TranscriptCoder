package au.org.intersect.samifier.domain;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import au.org.intersect.samifier.parser.TranscriptomeFileParsingException;

public class TranscriptInfo {

    private String transcriptId;
    private String chromosome;
    private String origin;          // Where the GTF file came from
    private int start;              // Start position of the transcript (WRT Genome)
    private int stop;               // End position of the transcript   (WRT Genome)
    private int direction;

    // Exon Table
    // ExonNumber|Chromosome|Strand|Start|End|Phase
    // 1         |1         |+     |     |   |
    // 2         |1         |+     |     |   |
    // 1         |2         |-     |     |   |
    // 2         |2         |-     |     |   |
    private Map<String, ExonInfo> exonMap;

    // Codon Table
    // ExonNumber|StartOrStop
    // 1         |Start
    // 1         |Stop
    // ----------------------
    // 1         |Start
    // 2         |Stop
    // ----------------------
    // 1         |Start
    // 2         |Start
    // 3         |Stop
    // 4         |Stop
    private Map<String, CodonInfo> codonMap;

    public TranscriptInfo(String transcriptId, String chromosome, String origin,
        int start, int stop, int direction) {

        this.transcriptId = transcriptId;
        this.chromosome   = chromosome;
        this.origin       = origin;
        this.start        = start;
        this.stop         = stop;
        this.direction    = direction;

        this.exonMap  = new HashMap<String, ExonInfo>();
        this.codonMap = new HashMap<String, CodonInfo>();
    }

    public TranscriptInfo(TranscriptInfo t) {
        this(t.getId(), t.getChromosome(), t.getOrigin(),
            t.getStart(), t.getStop(), t.getDirection());

        this.exonMap = new HashMap<String, ExonInfo>();
        for (ExonInfo oldE : t.getAllExons()) {
            ExonInfo newE = new ExonInfo(oldE);
            addExon(newE);
        }
    }

    // Get Functions
    public String getId() {
        return transcriptId;
    }
    
    public String getChromosome() {
        return chromosome;
    }

    public String getOrigin() {
        return origin;
    }

    public int getStart() {
        return start;
    }

    public int getStop() {
        return stop;
    }

    public int getDirection() {
        return direction;
    }

    public String getDirectionStr() {
        if (getDirection() == 1) {
            return GenomeConstant.FORWARD_FLAG;
        }
        return GenomeConstant.REVERSE_FLAG;
    }

    public Map<String, ExonInfo> getExonMap() {
        return exonMap;
    }

    public ExonInfo getExon(String exonId) {
        return exonMap.get(exonId);
    }

    public ExonInfo getExon(ExonInfo e) {
        List<ExonInfo> exons = getAllExons();
        int idx       = exons.indexOf(e);
        ExonInfo exon = exons.get(idx);
        return exon;
    }

    public List<ExonInfo> getAllExons() {
        List<ExonInfo> exons = new ArrayList<ExonInfo>(exonMap.values());
        Collections.sort(exons, new ExonInfoComparator());
        return exons;
    }

    public List<CodonInfo> getAllCodons() {
        List<CodonInfo> codons = new ArrayList<CodonInfo>(codonMap.values());
        Collections.sort(codons, new CodonInfoComparator());
        return codons;
    }

    public CodonInfo getStartCodon() {
        List<CodonInfo> codons = getAllCodons();
        if (codons.isEmpty()) {
            return null;
        }
        return codons.get(0);
    }

    public CodonInfo getStopCodon() {
        List<CodonInfo> codons = getAllCodons();
        if (codons.isEmpty()) {
            return null;
        }
        return codons.get(codons.size() - 1);
    }

    public List<String> getExonSequences(StringBuilder genomeString) {
        List<String> nucleotideSequences = new ArrayList<String>();
        for (ExonInfo e : getAllExons()) {
            String currNucleotideSequence = genomeString.substring(e.getStart() - 1, e.getStop());
            if (!isForward()) {
                String invertedSequence = ProteinOutputter.invertNucleotideSequence(currNucleotideSequence.toString()); 
                currNucleotideSequence  = new StringBuilder(invertedSequence).reverse().toString();
            }
            nucleotideSequences.add(currNucleotideSequence);
        }
        return nucleotideSequences;
    }

    // Add Functions
    public void addExon(ExonInfo exon) {
        exonMap.put(exon.getId(), exon);
    }

    public void addCodon(CodonInfo codon) {
        String key = codon.getExonId() + ":" + codon.getType();
        codonMap.put(key, codon);
    }

    // Has Functions
    public boolean hasExon(ExonInfo exon) {
        if (exonMap.containsKey(exon.getId())) {
            return true;
        }
        return false;
    }

    public boolean hasCodon(CodonInfo codon) {
        String key = codon.getExonId() + ":" + codon.getType();
        if (codonMap.containsKey(key)) {
            return true;
        }
        return false;
    }

    public boolean containsStartCodon() {
        List<CodonInfo> codons = getAllCodons();
        for (CodonInfo c : codons) {
            if (c.isStart()) {
                return true;
            }
        }
        return false;
    }

    public boolean containsStopCodon() {
        List<CodonInfo> codons = getAllCodons();
        for (CodonInfo c : getAllCodons()) {
            if (c.isStop()) {
                return true;
            }
        }
        return false;
    }

    public boolean isForward() {
        return getDirection() == 1;
    }

    public boolean isValid() {
        // Transcripts should either:
        //     * Have Start and Stop codons (Reference)
        //     * No start and stop codons (Target)
        if ((containsStartCodon() && containsStopCodon()) ||
            (!containsStartCodon() && !containsStopCodon())) {
            return true;
        }
        return false;
    }

    // Set Functions
    public void setId(String transcriptId) {
        this.transcriptId = transcriptId;
    }

    // Misc Functions
    public void update()
            throws TranscriptomeFileParsingException {

        Iterator<ExonInfo> iter = getAllExons().iterator();
        CodonInfo startCodon    = getStartCodon();
        CodonInfo stopCodon     = getStopCodon();

        // Only need to update transcripts that have
        // Start and Stop codons
        if (startCodon != null && stopCodon != null) {
            ExonInfo prevExon = null;
            while (iter.hasNext()) {
                ExonInfo currExon = iter.next();
                float eId      = Float.parseFloat(currExon.getId());
                float startEId = Float.parseFloat(startCodon.getExonId()); 
                float stopEId  = Float.parseFloat(stopCodon.getExonId()); 

                /**********
                 * Assign a score to each exon
                 *     * -1 = Non-coding
                 *     *  1 = Coding
                 *     *  0 = Start/Stop Boundary
                 * Theory: Higher scores means the exon is more likely
                 *         to be translated and therefore more likely to be better
                 *         for inferring the start codon of novel transcripts.
                 **********/
                // Exon is a Start/Stop boundary
                if (eId == startEId || eId == stopEId) {
                    currExon.setScore(0);

                // Exon is completely within the coing region (i.e., coding) 
                } else if (eId > startEId && eId < stopEId) {
                    currExon.setScore(1);

                // Exon is completely outside the coding region (i.e., non-coding)
                } else if (eId < startEId || eId > stopEId) {
                    currExon.setScore(-1);

                } else {
                    String errorMessage = "Exon " + eId + " in transcript " + getId() + "has unexpected type.";
                    throw new TranscriptomeFileParsingException(errorMessage);
                }

                /**********
                 * Annotate each exon with:
                 *     * The coding frame of the stop position
                 *     * Start and/or Stop codon information
                 **********/
                /**********
                 * We assign coding frames to the following exons:
                 *     * Boundary (Start)
                 *     * Coding
                 * Boundary (Stop) exons do not neeed a phase
                 **********/

                // Exon contains both Start and Stop codons
                // Special cases whereby only one exon is translated
                if (eId == startEId && eId == stopEId) {
                    // Not sure if theres a problem here...
                    currExon.setCodingFrame(0);
                    currExon.setStartCodon(startCodon);
                    currExon.setStopCodon(stopCodon);

                // Exon contains only Start codon
                // Boundary (Start) exon
                } else if (eId == startEId && eId != stopEId) {
                    int codingFrame = getCodingFrame(currExon, startCodon);
                    currExon.setCodingFrame(codingFrame);
                    currExon.setStartCodon(startCodon);

                // Exon contains only Stop codon
                // Boundary (Stop) exon)
                } else if (eId != startEId && eId == stopEId) {
                    currExon.setStopCodon(stopCodon);

                // Exon does not contain Start and Stop codons
                // Coding / Non-coding exons
                } else if (eId != startEId && eId != stopEId) {
                    // Check whether the exon is Coding or Non-coding
                    if (currExon.getScore() == 1) {
                        int codingFrame = getCodingFrame(currExon, prevExon);
                        currExon.setCodingFrame(codingFrame);
                    }
                }
                prevExon = currExon;
            }
        }
    }

    public int getCodingFrame(ExonInfo currExon, CodonInfo startCodon) {
        // Do forward calculation by default
        int numBases = Math.abs((currExon.getStop() - startCodon.getStart() + 1));
        if (!isForward()) {
            numBases = Math.abs((startCodon.getStop() - currExon.getStart() + 1));
        }
        return numBases % GenomeConstant.BASES_PER_CODON;
    }

    public int getCodingFrame(ExonInfo currExon, ExonInfo prevExon) {
        // Number of bases in the previous exon and the
        // Number of bases in the current exon to make up a codon
        int peBases  = GenomeConstant.BASES_PER_CODON - prevExon.getCodingFrame();
        int ceBases  = peBases % GenomeConstant.BASES_PER_CODON;
        int numBases = Math.abs(((currExon.getStop() - currExon.getStart() + 1) - ceBases));
        return numBases % GenomeConstant.BASES_PER_CODON; 
    }

    @Override
    public String toString() {
        StringBuffer out = new StringBuffer();
        out.append(transcriptId + "\t" + chromosome + "\t" + origin);
        out.append("\t");
        out.append(getDirectionStr() + "\t" + start + "\t" + stop);
        out.append(System.getProperty("line.separator"));

        List<ExonInfo> exons = getAllExons();
        for (ExonInfo e : exons) {
            out.append("\t");
            out.append(e);
            out.append(System.getProperty("line.separator"));
        }

        List<CodonInfo> codons = getAllCodons();
        for (CodonInfo c : codons) {
            out.append("\t\t");
            out.append(c);                   
            out.append(System.getProperty("line.separator"));
        }

        return out.toString();
    }

}
