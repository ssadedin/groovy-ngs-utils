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

    @Test
    public void readGencode() {
        
        println "${new Date()} Reading gencode from $gencode...."
        Gencode gencode = new Gencode(gencode)
        gencode.progress = new ProgressCounter(withRate: true, withTime: true, timeInterval:5000)
        gencode.load()
        
        println "${new Date()} done"

        println "Read: $gencode"
        
        def dvl1 = gencode.getGeneRegion("DVL1")
        
        println "DVL1 location = $dvl1 with id $dvl1.id" 
        
        println "${dvl1.children.size()} Transcripts of DVL1: " 
        
        dvl1.children.each { Transcript tx ->
            println "Transcript: $tx"
        }
        
        assert dvl1.children.size() == 2
        
        println "Exons of the main transcript: "
        
        dvl1.children[0].children.each {
            println "Exon $it"
        }
        
        assert dvl1.children[0].children.size() == 17
        assert dvl1.children[1].children.size() == 17
    }
    
    @Test
    void 'read indexed file'() {
        Gencode gencode = new Gencode(gencode)
        
        gencode.loadRegion(new Region('chr1:1335276-1349418'))

        def dvl1 = gencode.genes['DVL1']
        
        println "DVL1 is at $dvl1"
        
        println "There are ${dvl1.children.size()} transcripts: " + dvl1.children
        
        println "The first transcript has " + dvl1.children[1].children.size() + " exons: " + dvl1.children[1].children.join('\n')
        
        // For debug
//        dvl1.children.eachWithIndex { tx, i ->
//            Regions dvl_regions = new Regions(tx.children)
//            dvl_regions.save("dvl_regions_${i}.bed", extra: { tx.id })
//        }

        assert dvl1.children[0].children.size() == 17
        assert dvl1.children[1].children.size() == 17
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
