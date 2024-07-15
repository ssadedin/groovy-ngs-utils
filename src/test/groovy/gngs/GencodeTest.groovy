package gngs

import static org.junit.Assert.*

import gngs.gencode.Exon
import gngs.gencode.Gencode
import gngs.gencode.Transcript
import htsjdk.tribble.gff.Gff3Codec
import htsjdk.tribble.readers.AsciiLineReader
import htsjdk.tribble.readers.LineIteratorImpl
import org.junit.Test

class GencodeTest {
    
    String gencode = 'src/test/data/gencode.test.gff3.gz'
    
    Region DVL1 = new Region('chr1:1335276-1349418')

    @Test
    public void readGencode() {
        
        println "${new Date()} Reading gencode from $gencode...."
        Gencode gencode = new Gencode(gencode)
        gencode.progress = new ProgressCounter(withRate: true, withTime: true, timeInterval:5000)
        gencode.load()
        
        println "${new Date()} done"

        println "Read: $gencode"
        
        def dvl1 = gencode.getGeneFeature("DVL1")
        
        println "DVL1 location = $dvl1 with id $dvl1.id" 
        
        println "${dvl1.children.size()} Transcripts of DVL1: " 
        
        dvl1.children.each { Transcript tx ->
            println "Transcript: $tx"
        }
        
        assert dvl1.children.size() == 2
        
        println "Exons of the main transcript: "
        
        dvl1.children[0].children.each {
            println "Exon $it ($it.exonNumber, coding=$it.coding)"
        }
        
        assert dvl1.strand == '-'
        assert dvl1.children[0].children.size() == 17
        assert dvl1.children[1].children.size() == 17
        
        // DVL a case where coding sequence ends within the exon,
        // hence exon 15 represented twice, once with coding sequence,
        // once without
        assert dvl1.children[0].children[0].exonNumber == 15
    }
    
    @Test
    void 'read indexed file'() {
        Gencode gencode = new Gencode(gencode)
        
        gencode.loadRegion(DVL1)

        def dvl1 = gencode.genes['DVL1']
        
        println "DVL1 is at $dvl1"
        
        println "There are ${dvl1.children.size()} transcripts: " + dvl1.children
        
        println "The first transcript has " + dvl1.children[1].children.size() + " exons: " + dvl1.children[1].children.join('\n')
        
        // For debug
//        dvl1.children.eachWithIndex { tx, i ->
//            Regions dvl_regions = new Regions(tx.children)
//            dvl_regions.save("dvl_regions_${i}.bed", extra: { tx.id })
//        }

        assert dvl1.strand == '-'
        assert dvl1.children[0].children.size() == 17
        assert dvl1.children[1].children.size() == 17
    }
    
    
    @Test
    void 'gene symbols are resolved correctly'() {
        Gencode gencode = new Gencode(gencode).load()
        
        List<String> overlappingGenes = gencode.getGenes(DVL1)

        println "Gene symbols overlapping $DVL1 are " + overlappingGenes
        
        assert ['TAS1R3', 'DVL1', 'MIR6808'].every { overlappingGenes.contains(it) }
    }
    
    @Test
    void 'flattened exons are resolved correctly'() {
        Gencode gencode = new Gencode(gencode).load()
        Regions exons = gencode.getExons('DVL1')
        assert exons.numberOfRanges == 15
    }
    
    @Test
    void 'cds calculated correctly'() {
        Gencode gencode = new Gencode(gencode).load()
        
        def cds = gencode.getCDS(DVL1)
        
        println "CDS for DVL1 is $cds"
        
        def expected = [TAS1R3:0, DVL1:2073, MIR6808:0]
        assert cds.size() == 3
        assert cds.DVL1 == expected.DVL1
        assert cds.TAS1R3 == 0
        assert cds.MIR6808 == 0
        
    }

//    Only works on the full gencode
//
//    @Test
//    void 'read scn5a'() {
//        Gencode gencode = new Gencode(gencode)
//        
//        gencode.loadRegion(new Region('chr3:38,546,062-38,651,687'))
//        
//        def gene_symbol = 'SCN5A'
//
//        def scn5a = gencode.genes[gene_symbol]
//        
//        println "scn5a is at $scn5a"
//        
//        println "The transcripts are " + scn5a.children
//        
//        println "The first transcript has " + scn5a.children[1].children.size() + " exons: " + scn5a.children[1].children.join('\n')
//        
//        
////        scn5a.children.eachWithIndex { tx, i ->
////            Regions gene_regions = new Regions(tx.children)
////            gene_regions.save("${gene_symbol}_regions_${i}.bed", extra: { tx.id })
////        }
//    }
}
