import static org.junit.Assert.*;

import org.junit.Test;

class VariantDBTest {

    @Test
    public void testCountObservations() {
        def db = new VariantDB("small.db")
        
        def peds = Pedigrees.parse("/Users/simon/work/sjia/samples.ped")
        
        // find any old variant
        VCF vcf = VCF.parse("/Users/simon/work/sjia/variants/small.vcf", peds)
        
        Variant v = vcf[0]
        
        println "Observations for $v are: " + db.countObservations(v, v.alleles[0])
        
        println "Variant $v is in " + vcf.samples.count { v.sampleDosage(it) } + " samples"
        println "Variant $v is in " + v.pedigrees.size()  + " families"
    }

}
