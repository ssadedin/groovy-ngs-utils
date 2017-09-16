import java.nio.file.Files;

import com.xlson.groovycsv.PropertyMapper;

import groovy.transform.CompileStatic
import groovy.util.logging.Log;

/**
 * A convenience class that is set up to parse the UCSC RefGene database
 * as a RangedData object. This allows querying it over ranges to find 
 * transcripts and easy access to the properties of ranges so found.
 * <p>
 * A set of utility operations are additionally available for querying 
 * other features (lookup by gene, transcript name, etc).
 * 
 * @author simon
 */
@Log
class RefGenes {
    
    RangedData refData = null
    
    /**
     * Transcripts indexed by gene
     */
    Map<String,List<Integer>> geneToTranscripts = new HashMap(30000)
    
    Map<String, Region> transcriptIndex = new HashMap(80000)
    
    static COLUMN_NAMES= [ "num", "tx", "chr", "strand", "tx_start", "tx_end", "cds_start", "cds_end", "exons", "starts", "ends", 
                           "u1","gene", "cdsStartStat","cdsEndStat","exonFrames"]
    
    RefGenes(Map options = [:], String sourceFile) {
       this.refData = new RangedData(sourceFile, 2,4,5)
       this.refData.load(options + [columnNames:COLUMN_NAMES, zeroBased:true, readFirstLine:true])
       this.index() 
    }
    
    RefGenes(Reader r) {
        load(r)
    }
    
    /**
     * Default UCSC download URL for the refgene database
     */
    static String UCSC_REFGENE_URL = "http://hgdownload.soe.ucsc.edu/goldenPath/##genomeVersion##/database/refGene.txt.gz" 
    
    static RefGenes download(String genomeVersion="hg19") {
        
        Map<String,String> genomeMap = [
            "GRCh37" : "hg19",
            "GRCh38" : "hg38"
        ]
        
        boolean stripChr = false
        
        // Map to appropriate UCSC genome and strip chr if necessary
        String ucscGenomeVersion = genomeVersion
        if(genomeVersion in genomeMap) {
            ucscGenomeVersion = genomeMap[genomeVersion]
            stripChr = true
        }
        
        File outputFile = new File("refGene.txt.gz")
        if(!outputFile.exists()) {
            
            String ucscUrl = UCSC_REFGENE_URL.replace('##genomeVersion##',ucscGenomeVersion)
            outputFile.withOutputStream { outputStream ->
                new URL(ucscUrl).withInputStream { urlStream ->
                    Files.copy(urlStream, outputStream)
                }
            }
        }
        return new RefGenes(outputFile.name, stripChr:stripChr)
    }
    
    void load(Reader r) {
       this.refData = new RangedData(r, 2,4,5)
       this.refData.load(columnNames:COLUMN_NAMES, zeroBased:true, readFirstLine:true)
       this.index()
    }
    
    @CompileStatic
    void index() {
       int i=0
       for(refLine in this.refData) {
           String gene = (String)refLine.getProperty('gene')
           List<Integer> txes = (List<Integer>)geneToTranscripts[gene]    
           if(!txes)
               geneToTranscripts[gene] = txes = []
           txes.add(i)    
           ++i
           transcriptIndex[(String)refLine.getProperty('tx')] = refLine
       }
    }
    
    /**
     * Return the names of all the transcripts (NM_... identifiers) for the
     * given HGNC symbol
     * 
     * @param gene
     * @return
     */
    List<String> getTranscriptIds(String gene) {
        geneToTranscripts[gene].collect { refData[it].tx }
    }
    
    /**
     * Return a list of transcript objects as Regions having all the 
     * full properties of a transcript (see COLUMN_NAMES)
     */
    List<Region> getTranscripts(String gene) {
        geneToTranscripts[gene].collect { refData[it] }
    }
    
    /**
     * Return a regions object containing a region for each unique exon in the transcript.
     * If two exons overlap, the union of the two is returned.
     * 
     * @param gene
     * @return
     */
    @CompileStatic
    Regions getExons(String gene, boolean codingOnly=true) {
        Regions exons = new Regions()
        String strand
        geneToTranscripts[gene]?.collect { (Region)refData[it] }?.each { Region tx ->
            
            // Transcripts starting with NR are non-coding
            if(codingOnly && ((String)tx['tx']).startsWith("NR_"))
                return null
            
            strand = tx.getProperty('strand')
                
            [ ((String)tx['starts']).tokenize(","), ((String)tx['ends']).tokenize(",") ].transpose().collect { se ->
                
                List<String> exonStartEnd = (List<String>)(se);
                int start = exonStartEnd[0].toInteger()
                int end = exonStartEnd[1].toInteger()
                if(codingOnly) {
                    start = Math.max(((String)tx['cds_start']).toInteger(), start)
                    end = Math.min(((String)tx['cds_end']).toInteger(), end)
                }
                if(start>end)
                    (Region)null
                else {
                    new Region(tx.chr, new GRange(start+1,end+1, null))
                }
            }.each { Region r ->
                if(r != null)
                    exons.addRegion(r)
            }
        }
        
        Regions results = exons.reduce().enhance()
        results.each { it.setProperty('strand', strand)}
        return results
    }
    
    Regions getTranscriptExons(String transcriptId) {
        Region tx = this.transcriptIndex[transcriptId]
        def result = new Regions()
        [ tx.starts.split(","), tx.ends.split(",") ].transpose().each { 
           result.addRegion(new Region(tx.chr, (it[0].toInteger()+1)..(it[1].toInteger()+1)))
        }
        return result
    }
    
    /**
     * Return a regions object containing a region for each unique exon in the transcript.
     * If two exons overlap, the union of the two is returned.
     * 
     * @param gene
     * @return
     */
    Regions getSpliceSites(String transcript) {
        Regions splices = new Regions()
        def tx = transcriptIndex[transcript]
        if(tx==null)
            throw new IllegalArgumentException("No transcript $transcript could be found")
            
        [ tx.starts.split(","), tx.ends.split(",") ].transpose().collect { 
                new Region(tx.chr, (it[0].toInteger()+1)..(it[1].toInteger()+1)) 
        }.each {
            splices.addRegion(it)
        }
        splices[0].from = tx.cds_start+1
        splices[-1].to = tx.cds_end+1
        
        return splices
    }
    
    /**
     * Look up the entire span of a gene from its HGNC symbol
     * 
     * @param hgncSymbol
     * @return
     */
    Region getGeneRegion(String hgncSymbol) {
       Regions exons = getExons(hgncSymbol, false) 
       if(exons.numberOfRanges == 0)
           return null // unknown gene
       
       new Region(exons[0].chr, exons*.from.min()..exons*.to.max())
    }
    
    List<String> getGenes(Region r) {
        refData.getOverlaps(r)*.extra*.gene.unique()
    }
}