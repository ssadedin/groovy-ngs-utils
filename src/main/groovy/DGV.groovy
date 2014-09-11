
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
     
    DGV load() {
        this.dgv = Utils.time("Loading DGV data ...") { new RangedData(dgvFile, 1,2,3).load(columnNames:DGV_COLUMNS) }
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
}
