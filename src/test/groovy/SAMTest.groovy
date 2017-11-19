import static org.junit.Assert.*;

import org.junit.Test

import gngs.Regions
import gngs.SAM
import gngs.Utils
import groovy.transform.CompileStatic
import htsjdk.samtools.SAMRecord;;;


class SAMTest {

//    @Test
    public void testHardClipped() {
        
        SAM sam = new SAM("tests/data/hard_clipped.bam")
        
        sam.pileup("chr12", 53819703, 53819708) { PileupIterator.Pileup p ->
            println p.position + " : " + p.alignments.size()
        }
        
//        sam.coverage("chr9",98211125,98211141) {  p ->
//             println p.position   
//        }
//        
//        
//        println sam.meanCoverage("chr9", 98211125, 98211141)
        
        println sam.meanCoverage("chr12", 53819703, 53819708)
    }
    
    @Test
    public void testWholeChromosomeCoverage() {
        SAM sam = new SAM("/Users/simon/work/dsd/batch4/work/ZU1123_S11_L001_R1_001.fastq.trim.atrim.reorder.realign.recal.bam")
        
        Regions regions = new BED("/Users/simon/work/dsd/batch3/design/target_regions.bed").load()
        regions = regions.grep { it.chr == "chrY" } as Regions
        
        def stats = sam.coverageStatistics(regions)
        
        println "Mean coverage = $stats.mean reads counted = $stats.n"
    }

} 
