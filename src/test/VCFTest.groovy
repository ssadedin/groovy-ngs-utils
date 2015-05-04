import static org.junit.Assert.*;

import org.junit.Test;


class VCFTest {
    
    String testVCF = 
    
"""##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=2014-12-22
##source=IlluminaPlatinumGenomes, version: 7.0.0-production
##reference=hg19
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12877\tNA12878
chr1\t325075\t.\tG\tC\t.\tInconsistent\tMTD=bwa_gatk3;HD=12;ED=14\tGT:LOOHD:PF\t.:10:.\t.:10:.
chr1\t325127\t.\tT\tC\t.\tInconsistent\tMTD=bwa_gatk3;HD=22;ED=23\tGT:LOOHD:PF\t.:20:.\t0/1:22:T
chr1\t325155\t.\tC\tA\t.\tInconsistent\tMTD=bwa_gatk3;HD=20;ED=21\tGT:LOOHD:PF\t.:18:.\t1/1:20:T
chr1\t664608\t.\tG\tC\t.\tQCFail\tMTD=bwa_gatk3;HD=0;ED=5\tGT:LOOHD:PF\t0|0:0:.\t0|0:0:.
chr1\t665058\t.\tC\tG\t.\tInconsistent\tMTD=bwa_gatk3;HD=3;ED=4\tGT:LOOHD:PF\t0/0:2:.\t0/0:3:.
chr1\t762273\t.\tG\tA\t.\tPASS\tMTD=bwa_gatk3;HD=0;ED=5\tGT:LOOHD:PF\t1|0:0:T\t1|1:0:T
chr1\t762472\t.\tC\tT\t.\tPASS\tMTD=bwa_gatk3;HD=0;ED=5\tGT:LOOHD:PF\t0|1:0:T\t0|0:0:.
chr1\t762485\t.\tC\tA\t.\tPASS\tMTD=bwa_gatk3;HD=0;ED=5\tGT:LOOHD:PF\t0|1:0:T\t0|0:0:.
"""
    
    @Test
    public void testInfoMetaData() {
        String info = """##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence type as predicted by VEP. Format: Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|CLIN_SIG|CANONICAL|AFR_MAF|AMR_MAF|ASN_MAF|EUR_MAF|SYMBOL|SYMBOL_SOURCE|PUBMED|AA_MAF|EA_MAF|PolyPhen|HGVSc|HGVSp|ENSP|DISTANCE|SIFT|Condel">"""
        VCF vcf = new VCF()
        vcf.headerLines = [info]
        
        def csqInfo = vcf.getInfoMetaData("CSQ")
        println "csqInfo = $csqInfo"
        
        def vepColumns = vcf.getVepColumns()
        println "Vep columns are $vepColumns"
        
        assert "Feature" in vepColumns
    }

    @Test
    public void testIn() {
        VCF vcf = new VCF()
        Variant v = new Variant(chr:"chr1", pos:10, ref:"C",  alt:"A", alts:["A","G"])
        vcf.add(v) 
        
        assert v in vcf
        
        Variant v3 = new Variant(chr:"chr2", pos:10, ref:"C", alt:"A")
        assert !(v3 in vcf) : "Variant with different chromosome still in VCF!"
        
        Variant v4 = new Variant(chr:"chr1", pos:11, ref:"C", alt:"A")
        assert !(v4 in vcf) : "Variant with different position still in VCF!"
        
        Variant v5 = new Variant(chr:"chr1", pos:10, ref:"C", alt:"T")
        assert !(v5 in vcf) : "Variant with different allele still in VCF!"
        
        Variant v6 = new Variant(chr:"chr1", pos:10, ref:"C", alt:"T")
        assert !(v6 in vcf) : "Variant with different reference still in VCF! (nonsensical anyway?)"
        
        Variant v7 = new Variant(chr:"chr1", pos:10, ref:"C", alt:"G")
        assert v7 in vcf : "Variant with multiple alternate alleles did not recognize second alternate"
            
        for(Variant v2 in vcf) {
            assert v2 == v
        }
    }
    
    @Test
    void testSelectSample() {
        VCF vcf = VCF.parse(samples: ["NA12878"], new ByteArrayInputStream(testVCF.bytes) )
        
        assert vcf.samples.size() == 1
        
        assert vcf[0].line.indexOf("NA12877") == -1
        
        vcf = VCF.parse(samples: ["NA12878"], new ByteArrayInputStream(testVCF.bytes))
        
        assert !vcf.variantsAt("chr1", 762273).empty
    }

}
