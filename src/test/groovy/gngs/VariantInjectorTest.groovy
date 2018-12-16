package gngs

import static org.junit.Assert.*

import org.junit.Test

class VariantInjectorTest {
    
    FASTA reference = new FASTA("src/test/data/small.fasta")

    @Test
    public void testInjectSNV() {
        VariantInjector inj = new VariantInjector(reference, var(30, 'T','C'))
        inj.alleleFraction = 1.0
        
        FASTQRead r1 = new FASTQRead('TATTAACCACTCACGGGAGCTC')
        inj.inject(r1, r1)
        
        assert inj.total > 0
        assert r1.bases == 'TATTAACCACCCACGGGAGCTC'
        
    }
    
    @Test
    public void testInjectDel() {
        VariantInjector inj = new VariantInjector(reference, var(30, 'TC','T'))
        inj.alleleFraction = 1.0
        
        FASTQRead r1 = new FASTQRead('TATTAACCACTCACGGGAGCTC')
        inj.inject(r1, r1)
        
//        println inj.referenceBases
//        println Align.global(inj.referenceBases, inj.replacementBases).profile
//        
        assert inj.total > 0
        assert r1.bases == 'TATTAACCACTACGGGAGCTC'
        
    }
    
    @Test
    public void testInjectIns() {
        VariantInjector inj = new VariantInjector(reference, var(30, 'T','TG'))
        inj.alleleFraction = 1.0
        
        FASTQRead r1 = new FASTQRead('TATTAACCACTCACGGGAGCTCT')
        inj.inject(r1, r1)
        
//        print Align.global(inj.referenceBases, inj.replacementBases).profile
//        println inj.referenceBases
//        println r1.bases
        
        assert inj.total > 0
        assert r1.bases == 'TATTAACCACTGCACGGGAGCTC'
    } 
    
    @Test
    public void testInjectRepeatIns() {
        
        // A/AAGAT
        
        FASTA ref = new FASTA() {
            String basesAt(String chr, long start, long end) {
                
                return 'CTGTCCTGAAAGATAGATATTAATTTCA'
//                return 'TGAAATTAATATCTATCTTTCAGGACAG'
            }
        }
        
        VariantInjector inj = new VariantInjector(ref, var(117188810, 'A','AAGAT'))
        inj.alleleFraction = 1.0
        
//        print Align.global(inj.referenceBases, inj.replacementBases).profile
        
        assert inj.replacementBases == 'CTGTCCTGAAAAGATGATAGATATTAATTTCA'
        
//        String readBases = 'CAAGCATCTATTGAAAATATCTGACAAACTCATCTTTTATTTTTGATGTGTGTGTGTGTGTGTGTGTGTTTTTTTAACAGGGATTTGGGGAATTATTTGAGAAAGCAAAACAAAACACTAACAATAGAAAAACTTCTAATGGTGATGACAGCCTCTTCTTCAGTAATTGCTCACTTCTTGGTACTCCTGTCCTGAAAGATATTAATTTCAAGATAGAAAGAGGACAGTTGTTGGCGGTTGCTGGATCCACTGGAGCAGGCAAGGTAGTTCTTTTTG'
//        FASTQRead r1 = new FASTQRead(readBases)
//        FASTQRead r2 = new FASTQRead(FASTA.reverseComplement(r1.bases))
//        inj.inject(r1, r2)
//        
//        println inj.referenceBases
//        println r1.bases
//        
//        assert inj.total > 0
    } 
    
    @Test
    public void testDelIns() {
        
        FASTA ref = new FASTA() {
            String basesAt(String chr, long start, long end) {
//                return 'AGTGCAGAAAGAAGAAATTCAATCCT'
                
                return 'TAGTGCAGAAAGAAGAAATTCAATCCTAA'
            }
        }
  
        VariantInjector inj = new VariantInjector(ref, var(13, 'GAAATTCAATCCT','AGAAA'))
        inj.alleleFraction = 1.0
        
        //                                         GAAATTCAATCCT
        FASTQRead r1 = new FASTQRead("AGTGCAGAAAGAAGAAATTCAATCCTAACTGAGACCTTA")
        inj.inject(r1, r1)
        
        print Align.global(inj.referenceBases, inj.replacementBases).profile
//        println inj.referenceBases
//        println r1.bases
        
        println inj.replacementBases
        println '-------GAAAGAAAGAAAAA'
        //       GAAAGAAAGAAAAACT
        //       GAAAGAAAGAAAAA

        assert inj.replacementBases.indexOf('GAAAGAAAGAAAA') >= 0
    }
    
    
    FASTQRead fq(String bases) {
        new FASTQRead(bases)
    }
    

    
    Variant var(int pos, String ref, String alt) {
        VCF vcf = new VCF()
        vcf.lastHeaderLine = ['CHROM','POS','REF','ALT','QUAL','FILTER','INFO','FORMAT','JOHNSMITH'] as String[]
        vcf.headerLines = [
            '##fileformat=VCFv4.2',
            '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">',
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            vcf.lastHeaderLine.join('\t')
        ]
        vcf.samples = ['JOHNSMITH']
        
        String line = "chr1 $pos    .  $ref $alt  32852.19    PASS    AC=1,1;AF=0.500,0.500;AN=2;DB;DP=158;FS=0.000;HaplotypeScore=84.0101;MLEAC=1,1;MLEAF=0.500,0.500;MQ=59.01;MQ0=0;QD=207.93;RPA=9,7,8;RU=AAG;STR;set=Intersection GT:AD:DP:GQ:PL  1/2:0,12,125:155:99:6991,6090,6305,609,0,197" 
        Variant v = Variant.parse(line.replaceAll(' {1,}','\t'))
        v.header = vcf
        return v
    }

}
