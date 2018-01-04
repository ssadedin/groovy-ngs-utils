import static org.junit.Assert.*;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import gngs.VCF
import gngs.Variant

class VariantDBTest {

    VCF vcf
    
    VariantDB db
    
    @Before
    void setup() {
        if(new File("test.db").exists())
            new File("test.db").delete()
        
        vcf = VCF.parse(new ByteArrayInputStream(VCFTest.testVCF.bytes))
        db = new VariantDB("test.db")
        vcf.each { db.add("batch1", null, it, null, "testcohort") }
        println "Added $vcf.size variants to database for samples $vcf.samples"
    }
    
    @After
    void after() {
        try {
            db.close()
        }
        catch(Exception e) {
            println "Exception closing db: $e"
        }
    }
    
    void addBatch(int index) {
        VCF vcf2 = new VCF(vcf)
        vcf2.samples= ["NA12877-"+index.toString(), "NA12878-"+index.toString()]
        
        // Add the new samples
        vcf2.each { db.add("batch$index", null, it, null, "testcohort") }
    }
    
    
    @Test
    void testCountVariants() {
        
        Variant v1 = vcf[0]
        Variant v2 = vcf[1]
        
        def obs = db.countObservations(v1, v1.alleles[0])
        println "Observations of $v1 are: " + obs
        
        // The first variant is in the VCF but is NOT genotyped to exist in either sample
        assert obs.samples == 0
        assert obs.families == 0
        
        
        // The second variant is only in 1 sample
        obs = db.countObservations(v2, v2.alleles[0])
        assert obs.samples == 1
        assert obs.families == 1
        
        Variant v3 = vcf.variantsAt("chr1", 762273)[0]
        obs = db.countObservations(v3, v3.alleles[0])
        assert obs.samples == 2
        assert obs.families == 2 // 2 families because we did not supply PED
    }
    
    @Test
    void testIgnoreNewVariants() {
        
        Thread.sleep(2000) // enough to ensure later batch time
        
        addBatch(2)
        
        Variant v2 = vcf[1]
        def obs = db.countObservations(v2, v2.alleles[0])
        
        assert obs.samples == 2
        assert obs.families == 2
        
        // Do not count observations that were added after 
        // the batch that the given sample was added
        obs = db.countObservations(v2, v2.alleles[0], "batch1")
        
        assert obs.samples == 1
        assert obs.families == 1
     }
    
    @Test
    void testQueryVariantCounts() {
        
        Variant v1 = vcf[0]
        Variant v2 = vcf[1]
        Variant v3 = vcf.variantsAt("chr1", 762273)[0]
        
        assert v3 != null
         
        // v1 is was not genotyped in any sample
        def obs = db.queryVariantCounts(v1, v1.alleles[0], null, "testcohort")
        assert obs.total == 0
        assert obs.in_target == 0
        assert obs.other_target == 0
        
        // v2 genotyped as present in 1 sample
        obs = db.queryVariantCounts(v2, v2.alleles[0], null, "testcohort")
        assert obs.total == 1
        assert obs.other_target == 0
        assert obs.in_target == 1
        
        // different cohort should move it out of target
        obs = db.queryVariantCounts(v2, v2.alleles[0], null, "testcohort2")
        assert obs.total == 1
        assert obs.other_target == 1
        assert obs.in_target == 0
    }
    
    @Test
    void testQueryByBatch() {
        
        Variant v1 = vcf[0]
        Variant v2 = vcf[1]
        Variant v3 = vcf.variantsAt("chr1", 762273)[0]
        
        assert v3 != null
        
        Thread.sleep(2000)
        
        addBatch(2)
         
        // v2 genotyped as present in 2 samples
        def obs = db.queryVariantCounts(v2, v2.alleles[0], null, "testcohort")
        assert obs.total == 2
        assert obs.other_target == 0
        assert obs.in_target == 2
        
        // If we specify batch1 we shoud only get 1
        obs = db.queryVariantCounts(v2, v2.alleles[0], null, "testcohort", batch:'batch1')
        assert obs.total == 1
        assert obs.other_target == 0
        assert obs.in_target == 1
         
        // different cohort should move it out of target
        obs = db.queryVariantCounts(v2, v2.alleles[0], null, "testcohort2")
        assert obs.total == 2
        assert obs.other_target == 2
        assert obs.in_target == 0
        
        // specify batch 1, should only see 1
        obs = db.queryVariantCounts(v2, v2.alleles[0], null, "testcohort2", batch:'batch1')
        assert obs.total == 1
        assert obs.other_target == 1
        assert obs.in_target == 0 
        
        Thread.sleep(1000)
        addBatch(3)
        
        Thread.sleep(1000)
        
        // should still only see 1
        obs = db.queryVariantCounts(v2, v2.alleles[0], null, "testcohort2", batch:'batch1')
        assert obs.total == 1
        assert obs.other_target == 1
        assert obs.in_target == 0 
         
        // now batch 2 should see 2
        obs = db.queryVariantCounts(v2, v2.alleles[0], null, "testcohort2", batch:'batch2')
        assert obs.total == 2
        assert obs.other_target == 2
        assert obs.in_target == 0     
    
        // batch 3 should see 3
        obs = db.queryVariantCounts(v2, v2.alleles[0], null, "testcohort2", batch:'batch3')
        assert obs.total == 3
        assert obs.other_target == 3
        assert obs.in_target == 0     
      }
    
    @Test
    void testReAddDifferentCohort() {
        
        Variant v2 = vcf[1]
        
        // Query same cohort - standard case, should come back in target
        def obs = db.queryVariantCounts(v2, v2.alleles[0], null, "testcohort")
        assert obs.total == 1
        assert obs.other_target == 0 
        assert obs.in_target == 1 
        
        VCF vcf2 = new VCF(vcf)
        vcf2.samples= ["NA12877", "NA12878"] // Same sample id but different cohort
        vcf2.each { db.add("batch2", null, it, null, "differentcohort") }    
        
        // Query the same cohort again - this should still come back 'intarget'
        obs = db.queryVariantCounts(v2, v2.alleles[0], null, "testcohort")
        assert obs.total == 1
        assert obs.other_target == 0
        assert obs.in_target == 1 
        
    }
}
