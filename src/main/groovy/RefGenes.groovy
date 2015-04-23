import com.xlson.groovycsv.PropertyMapper;

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
    
    RefGenes(String sourceFile) {
       this.refData = new RangedData(sourceFile, 2,4,5)
       this.refData.load(columnNames:COLUMN_NAMES, zeroBased:true, readFirstLine:true)
       
       this.index()
    }
    
    void index() {
       log.info "Indexing genes ..."
       int i=0
       for(refLine in this.refData) {
           List txes = geneToTranscripts[refLine.gene]    
           if(!txes)
               geneToTranscripts[refLine.gene] = txes = []
           txes.add(i)    
           ++i
           transcriptIndex[refLine.tx] = refLine
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
    Regions getExons(String gene, codingOnly=true) {
        Regions exons = new Regions()
        geneToTranscripts[gene].collect { refData[it] }.collect { tx ->
            
            // Transcripts starting with NR are non-coding
            if(!codingOnly && tx.tx.startsWith("NR_"))
                return null
                
            [ tx.starts.split(","), tx.ends.split(",") ].transpose().collect { exonStartEnd ->
                int start = exonStartEnd[0].toInteger()
                int end = exonStartEnd[1].toInteger()
                if(!codingOnly) {
                    start = Math.max(tx.cds_start.toInteger(), start)
                    end = Math.min(tx.cds_end.toInteger(), end)
                }
                if(start>end)
                    null
                else
                    new Region(tx.chr, (start+1)..(end+1)) 
            }.each {
                if(it != null)
                    exons.addRegion(it)
            }
        }
        exons.reduce()
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
}