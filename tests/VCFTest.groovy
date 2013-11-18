import static org.junit.Assert.*;

import org.junit.Test;


class VCFTest {

    @Test
    public void testIn() {
        VCF vcf = new VCF()
        Variant v = new Variant(chr:"chr1", pos:10, alt:"A", ref:"C", alts:["A","G"])
        vcf.add(v) 
        
        assert v in vcf
        
        Variant v3 = new Variant(chr:"chr2", pos:10, alt:"A", ref:"C")
        assert !(v3 in vcf) : "Variant with different chromosome still in VCF!"
        
        Variant v4 = new Variant(chr:"chr1", pos:11, alt:"A", ref:"C")
        assert !(v4 in vcf) : "Variant with different position still in VCF!"
        
        Variant v5 = new Variant(chr:"chr1", pos:10, alt:"T", ref:"C")
        assert !(v5 in vcf) : "Variant with different allele still in VCF!"
        
        Variant v6 = new Variant(chr:"chr1", pos:10, alt:"T", ref:"G")
        assert !(v6 in vcf) : "Variant with different reference still in VCF! (nonsensical anyway?)"
        
        Variant v7 = new Variant(chr:"chr1", pos:10, alt:"G", ref:"C")
        assert v7 in vcf : "Variant with multiple alternate alleles did not recognize second alternate"
            
        for(Variant v2 in vcf) {
            assert v2 == v
        }
        
        
        
    }

}
