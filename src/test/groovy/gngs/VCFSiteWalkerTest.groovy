package gngs

import static org.junit.Assert.*

import org.junit.Test

class VCFSiteWalkerTest {
    
    // Test VCF assumptions
    // pos 16591593 in both giab1, giab2
    // pos 16592176 only in giab2
    // pos 17209519 only in giab1
    // pos 17818113 G/A in giab2, G/C in giab1
    // pos 17906782 twice in giab2 as G/C, G/T, once in giab1
    // giab2 has overlaping indel - pos 18124323,18124324 

    @Test
    public void testWalkBySite() {
        
        def observed = []

        new VCFSiteWalker(['src/test/data/giab1.tiny.trio.vcf','src/test/data/giab2.tiny.trio.vcf']).walk { variants ->
            def pos = variants[0].start
            observed << pos
            switch(pos) {
                
                // pos 16591593 in both giab1, giab2
                case 16591593:
                    assert variants.size() == 2
                    break
                    
                // pos 16592176 only in giab2
                case 16592176:
                    assert variants.size() == 1
                    break
                    
                // pos 17209519 only in giab1
                case 17209519:
                    assert variants.size() == 1
                    break
                    
                // pos 17818113 G/A in giab2, G/C in giab1
                case 16591593:
                    assert variants.size() == 2
                    assert variants[0].alleles[0].baseString == 'A'
                    assert variants[1].alleles[0].baseString == 'C'
                    break

                // pos 17906782 twice in giab2 as G/C, G/T, once in giab1
                case 17906782:
                    assert variants.size() == 3
                    break

                // giab2 has overlaping indel - pos 18124323,18124324, giab1 only has 18124323
                case 18124323:
                    assert variants.size() == 2
                    break
                case 18124324:
                    assert variants.size() == 1
                    break
            }
        }
    
        List tested_positions =  [16591593, 16592176, 17209519, 16591593, 17906782, 18124323, 18124324]

        tested_positions.each { assert(it in observed) }
        
        assert tested_positions.indexOf(18124323) < tested_positions.indexOf(18124324) 
    }

    
    @Test
    public void testWalkByAllele() {
        
        def observed = []

        new VCFSiteWalker(['src/test/data/giab1.tiny.trio.vcf','src/test/data/giab2.tiny.trio.vcf']).walkByAllele { variants ->
            def pos = variants[0].start
            observed << pos
            switch(pos) {
                
                // pos 16591593 in both giab1, giab2
                case 16591593:
                    assert variants.size() == 2
                    break
                    
                // pos 16592176 only in giab2
                case 16592176:
                    assert variants.size() == 1
                    break
                    
                // pos 17209519 only in giab1
                case 17209519:
                    assert variants.size() == 1
                    break
                    
                // pos 17818113 G/A in giab2, G/C in giab1
                case 16591593:
                    assert variants.size() == 1
                    assert variants[0].alleles[0].baseString in ['A','C'] 
                    break

                // pos 17906782 twice in giab2 as G/C, G/T, once in giab1
                case 17906782:
                    assert variants.size() == 2 || variants.size() == 1
                    break

                // giab2 has overlaping indel - pos 18124323,18124324, giab1 only has 18124323
                case 18124323:
                    assert variants.size() == 2
                    break
                case 18124324:
                    assert variants.size() == 1
                    break
            }
        }
    
        List tested_positions =  [16591593, 16592176, 17209519, 16591593, 17906782, 18124323, 18124324]

        tested_positions.each { assert(it in observed) }
        
        assert tested_positions.indexOf(18124323) < tested_positions.indexOf(18124324) 
    } 
    
    @Test
    void testNamedVCFs() {

        def observed = []

        new VCFSiteWalker(['src/test/data/giab1.tiny.trio.vcf','src/test/data/giab2.tiny.trio.vcf'], ['joe','frank']).walkByAllele { variants ->
            int pos = variants[0].start
            if(pos == 16591593) {
                assert variants[0].source == 'joe'
                assert variants[1].source == 'frank'
            }
        }
    }
}
