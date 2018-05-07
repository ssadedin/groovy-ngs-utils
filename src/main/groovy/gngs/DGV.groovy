import java.util.zip.GZIPInputStream

import gngs.RangedData
import gngs.Region
import gngs.Utils
import groovy.transform.CompileStatic

/**
 * A thin wrapper around a RangedData object to add some utility functions
 * for accessing DGV data.
 * 
 * @author simon
 */
class DGV {
    
    /**
     * Columns from schema of DGV in UCSC table
     */
    List DGV_COLUMNS = ["bin", "chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "varType", "reference", "pubMedId", "method", "platform", "mergedVariants", "supportingVariants", "sampleSize", "observedGains", "observedLosses", "cohortDescription", "genes", "samples"]
    
    String dgvFile 
    
    @Delegate
    RangedData dgv
    
    DGV(String dgvFile) {
        this.dgvFile = dgvFile
    } 
    
    DGV(RangedData dgv) {
        this.dgv = dgv
    } 
     
    DGV parse() {
        this.dgv = new RangedData(dgvFile, 1,2,3).load(columnNames:DGV_COLUMNS) 
        return this
    }
    
    /**
     * Convenience method to return CNVS overlapping the specified region
     * 
     * @param r
     * @return
     */
    List<Region> queryOverlapping(Region r) {
       return dgv.getOverlaps(r)*.extra
    }
    
    /**
     * Find the maximum frequency of this CNV within any study within DGV where the 
     * study has more than a minimum threshold size (default: 10 people).
     * 
     * @param region    region to search
     * @return  maximum frequency, or zero if no CNVs found
     */
    double maxFreq(Map options=[:], Region region) {
        int minSampleSize = options.minSampleSize?:10
        List<Region> overlappingEntries = this.queryOverlapping(region)
        
        
        Region maxEntry = overlappingEntries.grep {
            (it.sampleSize > minSampleSize) && (it.observedGains + it.observedLosses < it.sampleSize)
        }.max {
            (it.observedGains + it.observedLosses) / it.sampleSize
        }
        
        if(maxEntry == null)
            return 0.0d
            
        return (maxEntry.observedGains + maxEntry.observedLosses) / maxEntry.sampleSize
    }
}
