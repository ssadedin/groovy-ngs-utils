package gngs.coverage

import gngs.Region
import gngs.Regions
import gngs.SAM
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SAMRecordIterator
import htsjdk.samtools.SamReader

/**
 * Optimised class to count reads over a given set of target regions
 */
@Log
class TargetedCoverageCalculator {
    
    SAM  bam
    
    Regions regions

    public TargetedCoverageCalculator(SAM bam, Regions regions) {
        super();
        this.bam = bam;
        this.regions = regions;
    }
    
    /**
     * Return an integer array containing the number of reads overlapping each of the intervals
     * in the bam file.
     * @return
     */
    @CompileStatic
    int [] countReads() {

        Map<String,Integer> chrCounts = regions.countBy { it.chr }

        List<String> chrs = chrCounts*.key
        
        log.info "Processing chromosomes $chrs"
        log.info "Read counts are $chrCounts"
        
        int [] results = new int[(int)chrCounts*.value.sum()]
        
        Region currentRegion = null
        int currentRegionIndex = 0
        int currentChromosomeIndex = -1
        String currentChromosome = null
        
       
        for(chr in chrs) {
            int startPos = (int)(chrCounts.takeWhile { 
                chr != it.key 
            }*.value.sum()?:0)

            calculateChrCoverage(chr, results, startPos)
        }
        return results
    }
    
    /**
     * Count the number of reads overlapping each target region from the given chromosome
     * 
     * @param chr       chromosome
     * @param results   results array, will be populated from startPos
     * @param startPos
     */
    @CompileStatic
    void calculateChrCoverage(String chr, int [] results, int startRegionIndex) {
        Regions chrRegions = regions.grep { ((Region)it).chr == chr } as Regions
        SamReader reader = bam.newReader(fast:true)
        int regionIndex = startRegionIndex
        Iterator<Region> regionIterator = chrRegions.iterator()
        Region currentRegion = regionIterator.next()
        int currentRegionStart = currentRegion.from
        int currentRegionEnd = currentRegion.to
        try {
            SAMRecordIterator iter = reader.query(chr, 0, 0, false)
            while(iter.hasNext()) {
                SAMRecord r = iter.next()
                final int start = r.alignmentStart
                final int end = r.alignmentEnd
                
                if(start > currentRegionEnd) {
                    while(true) {
                        if(!regionIterator.hasNext())
                            break finished

                        currentRegion = regionIterator.next()
                        ++regionIndex
                        if(currentRegion.to > start) {
                            currentRegionEnd = currentRegion.to
                            currentRegionStart = currentRegion.from
                            break
                        }
                    }
                }
                
                // Start is before the end of current region
                // 
                // 1. region is ahead of read
                //
                // |--read---|
                //                   |---region--|
                // 
                // 2. region is overlapping read
                // |--read---|
                //       |---region--|
                if(end > currentRegionStart) {
                    ++results[regionIndex]
                }
            }
            
            finished: 
               return
        }
        finally {
            reader.close()
        }
    }
}
