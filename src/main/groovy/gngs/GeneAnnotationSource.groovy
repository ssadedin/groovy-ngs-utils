package gngs

/**
 * Abstract interface for gene annotation sources
 * <p>
 * Enables dependent code to work with different sources of 
 * gene annotations, eg: Ensembl vs RefSeq
 * 
 * @author simon.sadedin
 */
interface GeneAnnotationSource {
    
    /**
     * @return the list of all gene symbols overlapping the given region
     */
    List<String> getGenes(IRegion region) 
    
    /**
     * @return list of regions representing flattened exons (all transcripts) of the given gene
     */
    Regions getExons(String gene)
    
    /**
     * Compute a map of gene symbols with the amount of coding sequence overlapped by each gene 
     * for the given region
     * 
     * @return  Map from gene symbol to int representing size of CDS for that gene
     */
    Map<String,Integer> getCDS(final IRegion region) 
}
