package gngs

import static org.junit.Assert.*

import org.junit.Before
import org.junit.Test

class GnomADCNVDatabaseTest {
    
    Region duplicationRegion = new Region('chr1:11000-51000')
    
    Region deletionRegion = new Region('chr1:22000-30000')
    
    Region bndRegion = new Region('chr1:219643057-219644757')
    
    GnomADCNVDatabase gmd 
        
    @Before
    void loadCNVs() {
        gmd = new GnomADCNVDatabase(new VCFIndex('src/test/data/gnomad_sv_test.vcf.gz'))
    }

    @Test
    public void testParseNoBreakEnds() {
        assert gmd.queryOverlapping(bndRegion).size() == 0
    }

    @Test
    void 'test deletions are found'() {
        def dels =  gmd.queryOverlapping(deletionRegion)
        
        println "Found ${dels.size()} dels"

        assert dels.size() == 1
        assert dels[0].qual == 999
    }
    
    @Test
    void 'test duplications are found'() {
        def dups =  gmd.queryOverlapping(duplicationRegion)
        assert dups.size() == 3
        assert dups.any{it.qual == 999.0}
    } 
    
    @Test
    void 'test pop max'() {
        def dupMaxFreq =  gmd.maxFreq(duplicationRegion)
        println dupMaxFreq
        
        assert Math.abs(dupMaxFreq - 0.19) < 0.02
    }
    
    @Test
    void 'test contigs prefixed with chr are found'() {
        Region dupChrRegion = new Region(duplicationRegion.chr, duplicationRegion.from, duplicationRegion.to)
        def dups =  gmd.queryOverlapping(dupChrRegion)
        assert dups.size() == 3
        println dups*.qual
        assert dups.any{it.qual == 999.0}
    } 
 
}
