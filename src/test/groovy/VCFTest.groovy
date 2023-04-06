import static org.junit.Assert.*;

import org.junit.Test;
import gngs.ContigHeader
import gngs.FormatMetaData
import gngs.VCF
import gngs.Variant


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
    void testFindWithDifferentCase() {
        VCF vcf = new VCF()
        Variant v1 = new Variant(chr:"chr1", pos:10, ref:"C",  alt:"A", alts:["A","G"])
        vcf.add(v1) 
        assert v1 in vcf

        Variant v2 = new Variant(chr:"chr1", pos:10, ref:"c",  alt:"A", alts:["A","G"])
        assert v2 in vcf
        assert vcf.find(v2) != null
        
        Variant v3 = new Variant(chr:"chr1", pos:10, ref:"C",  alt:"a", alts:["a","G"])
        assert v3 in vcf
        assert vcf.find(v3) != null
        
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
        assert v in vcf
        
        v = new Variant(chr:"chr1", pos:1276973, ref:"G", alt:"GACAC")
        assert vcf.find(v) != null
        assert vcf.find(v).type == "INS"
        assert v in vcf
        
        v = new Variant(chr:"chr1", pos:38098177, ref:"CA", alt:"C")
        assert vcf.find(v) != null
        
        assert v in vcf
        
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
    
    @Test
    void testMergeBug() {
        
        def v1Text = 
"""##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=2014-12-22
##source=IlluminaPlatinumGenomes, version: 7.0.0-production
##reference=hg19
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12877
chr1\t911916\trs74045046\tC\tT\t2668.77\t.\tAC=1;AF=0.500;AN=2;BaseQRankSum=-0.721;DB;DP=239;FS=2.199;MLEAC=1;MLEAF=0.500;MQ=59.97;MQ0=0;MQRankSum=0.180;QD=11.17;ReadPosRankSum=-0.379;SOR=0.848;CSQ=T|ENSG00000187642|ENST00000341290|Transcript|synonymous_variant|1932|1896|632|R|agG/agA|rs74045046||||||C1orf170|HGNC|||ENSP00000343864|ENST00000341290.2:c.1896G>A|ENST00000341290.2:c.1896G>A(p.%3D)|0.26|0.14|0.07|0.16||,T|ENSG00000187583|ENST00000491024|Transcript|downstream_gene_variant||||||rs74045046|||671|||PLEKHN1|HGNC|||ENSP00000462558|||0.26|0.14|0.07|0.16||\tGT:AD:GQ:PL\t0/1:105,133:99:2697,0,1980
""".trim()

        def v2Text = """##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=2014-12-22
##source=IlluminaPlatinumGenomes, version: 7.0.0-production
##reference=hg19
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12877
chr6\t42626564\trs35713624\tT\tTA\t12.06\t.\tAC=1;AF=0.500;AN=2;BaseQRankSum=0.358;DB;DP=7;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=60.00;MQ0=0;MQRankSum=0.358;QD=1.72;ReadPosRankSum=-1.231;SOR=0.223;CSQ=A|ENSG00000024048|ENST00000372901|Transcript|splice_region_variant&intron_variant&feature_elongation||||||rs35713624||||||UBR2|HGNC|||ENSP00000361992|ENST00000372901.1:c.3242+2_3242+3insA|||||||\tGT:AD:GQ:PL\t0/1:2,3:19:49,0,19
""".trim()

      VCF v1 = VCF.parse(new ByteArrayInputStream(v1Text.bytes) )
      VCF v2 = VCF.parse(new ByteArrayInputStream(v2Text.bytes) )
      VCF result = v1.merge(v2)
      assert result.count { 1 } == 2
    }
    
    // @Test
    void testMergeBug2() {
        
        def dir = "/Users/simon/work/vcgs/verify_cluster_migration/merge_debug"
        
        VCF vcf1 = VCF.parse("$dir/1.tiny.vcf")
        VCF vcf2 = VCF.parse("$dir/2.tiny.vcf")
        
        VCF merged = vcf1.merge(vcf2)
        
        assert merged.variantsAt("chr18", 48137).size() == 1
    }
    
    @Test
    void testParseFormatLine() { 
        
        VCF vcf = new VCF()
        
        
        FormatMetaData meta = vcf.parseFormatMetaDataLine('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        
        assert meta.type == String
        
        println "description = " + meta.description 
        
        assert meta.description == "Genotype"
        
        assert meta.id == "GT"
        
    }
    
    @Test
    void testMergeDifferentVEPs() {
        
        VCF vcf1 = VCF.parse('tests/data/test.21.vcf')
        VCF vcf2 = VCF.parse('tests/data/test.23.vcf')
        
        VCF merged = vcf1.merge(vcf2)
        
        println "Merged VCF has ${merged.size()} variants"
        
        Variant v = merged[1]
        println "Merged vep info = " + v.info
        
        println "EA_MAF = " + v.vepInfo.EA_MAF
        println "EUR_MAF = " + v.vepInfo[0].EUR_MAF
        
        assert v.vepInfo[0].EUR_MAF.tokenize(':')[1].toFloat() > 0.9
        
        VCF merged2 = vcf2.merge(vcf1)
        
        v = merged[1]
        assert v.vepInfo[0].EUR_MAF.tokenize(':')[1].toFloat() > 0.9
        
    }
    
    @Test
    void 'empty VCF still has samples'() {
        VCF vcf = VCF.parse("src/test/data/empty_vcf.vcf")
        
        assert vcf.samples != null
        assert vcf.samples.size() > 0
        assert vcf.samples[0] != null
    }
    
    @Test
    void 'parse VCFs containing SVs'() {
        VCF vcf = VCF.parse("src/test/data/sv.vcf")
        assert vcf[0].size() > 0
        println vcf[0].type
    }
    
    @Test
    void 'parse VCFs containing SVs 2'() {
        VCF vcf = VCF.parse("/Users/simon.sadedin/work/ximmer/src/test/data/test.canvas.vcf")
        
        println vcf[0].size()
        assert vcf[0].size() > 100
        println vcf[0].type
    }
    
    @Test
    void 'test filtering with a discard writer'() {
        String includePath = 'src/test/data/includes.vcf'
        def discardPath = 'src/test/data/discards.vcf'
        new File(includePath).withOutputStream { o ->
            VCF.filter(new ByteArrayInputStream(testVCF.bytes), discardWriter: discardPath, outputStream: o) { Variant v ->
                v.pos < 665058
            }
        }
        
        assert new File(includePath).exists()
        assert new File(discardPath).exists()
        
        VCF incVcf = VCF.parse(includePath)
        VCF excVcf = VCF.parse(discardPath)
        
        assert incVcf.find { it.pos == 325075 }
        assert !incVcf.find { it.pos == 762485 }
        
        assert !excVcf.find { it.pos == 325075 }
        assert excVcf.find { it.pos == 762485 }
    }
    
    /*
    @Test
    void 'test filter header'() {
        VCF.filter('tests/data/test.23.vcf', updateHeader: { return it }) { 
            // noop
        }
    }
    
    @Test
    void 'test trioDenovoRate calculation'() {
        VCF trio = VCF.parse('test.trio.chr1.vcf')
        
        println trio.trioDenovoRate()
    }
    */

    
    @Test
    void 'parse VCF with LowDepth annotations'() {
        VCF vcf = VCF.parse('test.lowdepth.vcf')
    }

    @Test
    void 'parse VCF from File'() {
        VCF vcf = VCF.parse(new File('test.lowdepth.vcf')) {
            true
        }
        
        assert vcf
    }
    
    @Test
    void 'filter VCF from File'() {
        VCF.filter(new File('test.lowdepth.vcf')) {  
            true
        } 
        
        assert true
    }
    
    @Test
    void 'test filter from vcf file name'() {
        VCF.filter('src/test/data/tmp.mini.vcf') { false }
    }
    
    @Test
    void 'test filter from vcf file name gzipped'() {
        VCF.filter('test.vcf.gz') { false }
    }
	
	@Test
	void 'test parse contig header'() {
        String header = '##contig=<ID=chr17_gl000203_random,length=37498,assembly=hg19>'
		
        VCF vcf = VCF.parse('test.lowdepth.vcf')
		
		ContigHeader headerInfo = vcf.parseContigHeader(header)
		
		assert headerInfo .assembly == 'hg19'
		assert headerInfo.chr == 'chr17_gl000203_random'
		assert headerInfo.length == 37498
		
	}
 }
