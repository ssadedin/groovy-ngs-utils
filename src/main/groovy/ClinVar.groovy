import java.util.List;
import java.util.zip.GZIPInputStream

import gngs.RangedData
import gngs.Region
import gngs.Utils
import groovy.lang.Delegate;

/**
 * A thin wrapper around a RangedData object to add some utility functions
 * for accessing clinvar data.
 * 
 * @author simon
 */
class ClinVar {
   
    /**
     * Columns from schema of clinvar in UCSC table
     */
    List CLINVAR_COLUMNS = ["chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","reserved","blockCount","blockSizes","chromStarts","origName","type","geneId","geneSym","clinSign","snpId","nsvId","rcvAcc","testedInGtr","phenotype","origin","assembly","cytogenetic","reviewStatus","hgvsCod","hgvsProt","numSubmit","lastEval","guidelines","otherIds"]
    
    String clinvarFile 
    
    @Delegate
    RangedData clinvar
    
    ClinVar(String clinvarFile) {
        this.clinvarFile = clinvarFile
    } 
    
    ClinVar(RangedData clinvar) {
        this.clinvar = clinvar
    } 
     
    ClinVar parse() {
        this.clinvar = Utils.time("Loading ClinVar data ...") { 
            new RangedData(clinvarFile, 0,1,2).load(columnNames:CLINVAR_COLUMNS) 
        }
        return this
    }
    
    /**
     * Convenience method to return CNVS overlapping the specified region
     * 
     * @param r
     * @return
     */
    List<Region> queryOverlapping(Region r) {
       return clinvar.getOverlaps(r)*.extra
    }
}
