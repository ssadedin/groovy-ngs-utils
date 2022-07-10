import static org.junit.Assert.*;

import java.nio.file.FileSystem
import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Path

import org.junit.Test

import gngs.BED
import gngs.Regions
import gngs.SAM
import gngs.Utils
import gngs.PileupIterator
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
    
//    @Test
    public void testWholeChromosomeCoverage() {
        SAM sam = new SAM("/Users/simon/work/dsd/batch4/work/ZU1123_S11_L001_R1_001.fastq.trim.atrim.reorder.realign.recal.bam")
        
        Regions regions = new BED("/Users/simon/work/dsd/batch3/design/target_regions.bed").load()
        regions = regions.grep { it.chr == "chrY" } as Regions
        
        def stats = sam.coverageStatistics(regions)
        
        println "Mean coverage = $stats.mean reads counted = $stats.n"
    }
    
    @Test
    void testReadNioPath() {
        
//        File zipFile = new File('src/test/data/zipped_bams.zip')
//        URI uri = URI.create("jar:file:$zipFile.absolutePath")
//        
//        FileSystem fs = FileSystems.newFileSystem(uri, [create:"true"])
//        Path p = fs.getPath('foo/small.test.bam')
//        assert Files.exists(p)
//        
        SAM bam = new SAM(new File('src/test/data/small.test.bam').toPath())
        
        println "The contigs in the BAM are $bam.contigList"
        
        int i =0
        bam.eachRecord {  r -> 
            ++i
        } 
        
        println "There are $i reads in the BAM file"
        
        assert i == 1899
    }
    
    @Test
    void testWriteNioPath() {
        SAM bam = new SAM(new File('src/test/data/small.test.bam').toPath())
        def reads = []
        bam.eachRecord {  r -> 
            reads << r
        } 
         
        bam.withWriter('src/test/data/test.output.bam', false) { w ->
            reads.each {
                w.addAlignment(it)
            }
        }
        
        println "Wrote ${reads.size()} reads"
    }

} 
