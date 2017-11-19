import static org.junit.Assert.*;


import org.junit.Test;

import gngs.Regions
import gngs.SAM


class SexKaryotyperTest {

    @Test
    public void testKarotype() {
        SAM sam = new SAM("/Users/simon/work/dsd/batch4/work/ZU1123_S11_L001_R1_001.fastq.trim.atrim.reorder.realign.recal.bam")
        Regions regions = new BED("/Users/simon/work/dsd/batch3/design/target_regions.bed").load()
        
        println "All regions are $regions"
        
        def karyotyper = new SexKaryotyper(sam, regions)
        
        karyotyper.run()
        
        println "=" * 80
        println "xCoverage stats = " + karyotyper.xCoverage
        println "=" * 80
        println "yCoverage stats = " + karyotyper.yCoverage
        println "=" * 80
        println "autoCoverage stats = " + karyotyper.autosomeCoverage
        println "=" * 80
        println "Karyotyping results: " + karyotyper.sex

    }

}
