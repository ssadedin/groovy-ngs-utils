import static org.junit.Assert.*;

import org.junit.Test;


class VCFTest {
    
    static String testVCF = 
    
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
chr1\t1276973\trs70949572\tG\tGACAC\t6420.73\t.\tAC=2;AF=1.00;AN=2;DB;DP=95;FS=0.000;HaplotypeScore=20.2054;MLEAC=2;MLEAF=1.00;MQ=59.76;MQ0=0;QD=67.59;RPA=2,4;RU=AC;STR;EFF=INTRON(MODIFIER||||670|DVL1||CODING|NM_004421|5)\tGT:AD:DP:GQ:PL\t1/1:0,95:95:99:6458,286,0\t1/1:0,95:95:99:6458,286,0
chr1\t22446108\trs11397924\tC\tCT\t865.73\t.\tAC=1;AF=0.500;AN=2;BaseQRankSum=0.247;DB;DP=101;FS=4.039;HaplotypeScore=42.3410;MLEAC=1;MLEAF=0.500;MQ=59.77;MQ0=0;MQRankSum=0.070;QD=8.57;RPA=6,7;RU=T;ReadPosRankSum=-2.530;STR;EFF=UTR_3_PRIME(MODIFIER||||351|WNT4||CODING|NM_030761|5)\tGT:AD:DP:GQ:PL\t0/1:54,42:101:99:903,0,1294\t0/1:54,42:101:99:903,0,1294
chr1\t38098177\trs5773597\tCA\tC\t1179.73\t.\tAC=2;AF=1.00;AN=2;DB;DP=27;FS=0.000;HaplotypeScore=122.2078;MLEAC=2;MLEAF=1.00;MQ=59.14;MQ0=0;QD=43.69;EFF=INTRON(MODIFIER||||200|RSPO1||CODING|NM_001242910|1),INTRON(MODIFIER||||236|RSPO1||CODING|NM_001242909|1),INTRON(MODIFIER||||263|RSPO1||CODING|NM_001038633|1),INTRON(MODIFIER||||263|RSPO1||CODING|NM_001242908|1)\tGT:AD:DP:GQ:PL\t1/1:0,27:27:81:1217,81,0\t1/1:0,27:27:81:1217,81,0
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

    /*
chr1\t325075\t.\tG\tC\t.\tInconsistent\tMTD=bwa_gatk3;HD=12;ED=14\tGT:LOOHD:PF\t.:10:.\t.:10:.
chr1\t325127\t.\tT\tC\t.\tInconsistent\tMTD=bwa_gatk3;HD=22;ED=23\tGT:LOOHD:PF\t.:20:.\t0/1:22:T
chr1\t325155\t.\tC\tA\t.\tInconsistent\tMTD=bwa_gatk3;HD=20;ED=21\tGT:LOOHD:PF\t.:18:.\t1/1:20:T
chr1\t1276973\trs70949572\tG\tGACAC\t6420.73\t.\tAC=2;AF=1.00;AN=2;DB;DP=95;FS=0.000;HaplotypeScore=20.2054;MLEAC=2;MLEAF=1.00;MQ=59.76;MQ0=0;QD=67.59;RPA=2,4;RU=AC;STR;EFF=INTRON(MODIFIER||||670|DVL1||CODING|NM_004421|5)\tGT:AD:DP:GQ:PL\t1/1:0,95:95:99:6458,286,0\t1/1:0,95:95:99:6458,286,0
chr1\t22446108\trs11397924\tC\tCT\t865.73\t.\tAC=1;AF=0.500;AN=2;BaseQRankSum=0.247;DB;DP=101;FS=4.039;HaplotypeScore=42.3410;MLEAC=1;MLEAF=0.500;MQ=59.77;MQ0=0;MQRankSum=0.070;QD=8.57;RPA=6,7;RU=T;ReadPosRankSum=-2.530;STR;EFF=UTR_3_PRIME(MODIFIER||||351|WNT4||CODING|NM_030761|5)\tGT:AD:DP:GQ:PL\t0/1:54,42:101:99:903,0,1294\t0/1:54,42:101:99:903,0,1294
chr1\t38098177\trs5773597\tCA\tC\t1179.73\t.\tAC=2;AF=1.00;AN=2;DB;DP=27;FS=0.000;HaplotypeScore=122.2078;MLEAC=2;MLEAF=1.00;MQ=59.14;MQ0=0;QD=43.69;EFF=INTRON(MODIFIER||||200|RSPO1||CODING|NM_001242910|1),INTRON(MODIFIER||||236|RSPO1||CODING|NM_001242909|1),INTRON(MODIFIER||||263|RSPO1||CODING|NM_001038633|1),INTRON(MODIFIER||||263|RSPO1||CODING|NM_001242908|1)\tGT:AD:DP:GQ:PL\t1/1:0,27:27:81:1217,81,0\t1/1:0,27:27:81:1217,81,0
*/
    @Test
    void testFindVariant() {
        VCF vcf = VCF.parse(new ByteArrayInputStream(testVCF.bytes) )
        
        
        Variant v = new Variant(chr:"chr1", pos:325075, ref:"G", alt:"C")
        assert vcf.find(v) != null
        assert vcf.find(v).alt == "C"
        
        v = new Variant(chr:"chr1", pos:1276973, ref:"G", alt:"GACAC")
        assert vcf.find(v) != null
        assert vcf.find(v).type == "INS"
        
        v = new Variant(chr:"chr1", pos:38098177, ref:"CA", alt:"C")
        assert vcf.find(v) != null
        
    }
    
    
    @Test
    void testMergeVCFs() {
        
        String modifiedVCF = testVCF.replaceAll("NA12877","MA12877")
                                    .replaceAll("NA12878","MA12878")
                                    .split("\n")
                                    .grep { !it.contains("1276973") }
                                    .join("\n") + 
                                    "\nchr4\t55131001\trs61317836\tTA\tT\t120.73\t.\tAC=1;AF=0.500;AN=2;BaseQRankSum=-1.141;DB;DP=32;FS=0.000;HaplotypeScore=457.9188;MLEAC=1;MLEAF=0.500;MQ=34.44;MQ0=0;MQRankSum=1.902;QD=3.77;RPA=17,16;RU=A;ReadPosRankSum=-0.140;STR;EFF=INTRON(MODIFIER||||1089|PDGFRA||CODING|NM_006206|4)\tGT:AD:DP:GQ:PL\t0/1:18,13:32:99:158,0,325\t1/1:18,13:32:99:158,0,325" 
        
        VCF vcf1 = VCF.parse(new ByteArrayInputStream(testVCF.bytes) )
        VCF vcf2 = VCF.parse(new ByteArrayInputStream(modifiedVCF.bytes) )
        
        VCF merged = vcf1.merge(vcf2)
        
        Variant v = merged.find(new Variant(chr:"chr1", pos:1276973, ref:"G", alt:"GACAC"))
        
        assert v != null
        
        println v.line
        
        // Only in 1st VCF
        assert v.sampleDosage("NA12878") == 2
        assert v.sampleDosage("NA12877") == 2
        assert v.sampleDosage("MA12878") == 0
        assert v.sampleDosage("MA12877") == 0
        
        v = merged.find(new Variant(chr:"chr4", pos:55131001, ref:"TA", alt:"T"))
        
        // Only in 2nd VCF
        assert v != null
        assert v.sampleDosage("NA12877") == 0
        assert v.sampleDosage("NA12878") == 0
        assert v.sampleDosage("MA12877") == 1
        assert v.sampleDosage("MA12878") == 2
        
        v = merged.find(new Variant(chr:"chr1", pos:325127, ref:"T", alt:"C"))
        
        // In both
        assert v.sampleDosage("NA12878") == 1
        assert v.sampleDosage("MA12878") == 1
        
    }
 }
