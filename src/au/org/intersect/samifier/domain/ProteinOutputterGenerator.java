package au.org.intersect.samifier.domain;


public class ProteinOutputterGenerator implements
        ProteinLocationBasedOutputterGenerator {
    private String databaseName;
    private StringBuilder genomeString;
    private CodonTranslationTable translationTable;

    public ProteinOutputterGenerator(String databaseName,
            StringBuilder genomeString, CodonTranslationTable translationTable) {
        this.databaseName = databaseName;
        this.genomeString = genomeString;
        this.translationTable = translationTable;
    }
   
    public ProteinOutputter getOutputterFor(ProteinLocation proteinLocation) {
        // String databaseName, StringBuffer genomeString, CodonTranslationTable
        // translationTable
        return new ProteinOutputter(proteinLocation, databaseName,
                genomeString, translationTable);
    }
    
    // For generating the protein sequence of a transcript isoform
    public ProteinOutputter getOutputterFor(TranscriptInfo transcriptinfo) {
        // String databaseName, StringBuffer genomeString, CodonTranslationTable
        // translationTable
        return new ProteinOutputter(transcriptinfo, databaseName, 
        		genomeString, translationTable);
    }
}
