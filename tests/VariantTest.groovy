/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2013 Simon Sadedin, ssadedin<at>gmail.com
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

import static org.junit.Assert.*;

import org.junit.Test;


class VariantTest {

    Variant v
    String seq
    
    @Test
    public void testApplyDEL() {
        v = new Variant(ref:"CAG", alt:"C", pos: 20)
        
        seq = "TTTCCTTTTTCAGTGTCATGGAGGAACAGTGCTGTGAAAGGTGTCTCAGGGGCATGACTGCAAAGGCACAGCCCTTTCCTCCCTCTCTCTCGCAGGTGTGTG"
        
        def mutated = "TTTCCTTTTTCAGTGTCATGGAGGAACAGTGCTGTGAAAGGTGTCTCGGGCATGACTGCAAAGGCACAGCCCTTTCCTCCCTCTCTCTCGCAGGTGTGTG"
//        assert v.applyTo(seq, )
        
    }
    
// TODO: Fix this test, but needs an appropriate VCF header now
//    @Test
//    void testGetVepInfo() {
//        v = new Variant(ref:"CAG", alt:"C", pos: 20)
//        v.line = "chr1  57179502    .   G   T   51.8    PASS    AC=1;AF=0.5;AN=2;BaseQRankSum=-4.147;DP=42;Dels=0;FS=17.681;HaplotypeScore=72.6709;MLEAC=1;MLEAF=0.5;MQ=41.17;MQ0=0;MQRankSum=-2.887;QD=1.23;ReadPosRankSum=4.331;set=Intersection;CSQ=G|ENSG00000171456|ENST00000375687|Transcript|missense_variant|2937|2513|838|K/R|aAg/aGg|rs35632616||YES|0.03|0.01|0|0|ASXL1|HGNC||0.0335906|0.000116279|probably_damaging(0.977)|ENST00000375687.4:c.2513A>G|ENSP00000364839.4:p.Lys838Arg|ENSP00000364839||deleterious(0.05)|0.764   AD:DP:GQ:GT:PL  .:40:80:0/1:80,0,255"
//        println v.vepInfo
//    }
//    
    @Test void testFixGenotypeOrder() {
        v = new Variant(ref:"CAG", alt:"C", pos: 20)
        v.line = "chr1  57179502    .   G   T   51.8    PASS    AC=1;AF=0.5;AN=2;BaseQRankSum=-4.147;DP=42;Dels=0;FS=17.681;HaplotypeScore=72.6709;MLEAC=1;MLEAF=0.5;MQ=41.17;MQ0=0;MQRankSum=-2.887;QD=1.23;ReadPosRankSum=4.331;set=Intersection;EFF=DOWNSTREAM(MODIFIER||||728|C1orf168||CODING|NM_001004303|),UTR_3_PRIME(MODIFIER||||552|PRKAA2||CODING|NM_006252|9)   AD:DP:GQ:GT:PL  .:40:80:0/1:80,0,255"
        v.fixGenoTypeOrder()
        
        println v.line
        
        assert v.line.split('\t')[8] == 'GT:AD:DP:GQ:PL'
        assert v.line.split('\t')[9] == '0/1:.:40:80:80,0,255'
    }
    
    @Test
    void testAnnovarEqual() {
        
        v = Variant.parse("chr1 12908064    .   CA  C   2444.73 PASS    AC=1;AF=0.500;AN=2;BaseQRankSum=-2.337;DP=273;FS=1.058;HaplotypeScore=340.4396;MLEAC=1;MLEAF=0.500;MQ=46.84;MQ0=0;MQRankSum=2.472;QD=8.96;RPA=2,1;RU=A;ReadPosRankSum=-0.877;STR;set=Intersection   GT:AD:DP:GQ:PL  0/1:180,77:264:99:2482,0,6848")
        assert v.equalsAnnovar("chr1", 12908065, "-")
        
        v = Variant.parse("chr3 40503520    .   ACTGCTG ACTGCTGCTGCTGCTG,A  4028.19 PASS    AC=1,1;AF=0.500,0.500;AN=2;DP=70;FS=0.000;HaplotypeScore=190.6420;MLEAC=1,1;MLEAF=0.500,0.500;MQ=49.60;MQ0=0;QD=57.55;RPA=10,13,8;RU=CTG;STR;set=Intersection   GT:AD:DP:GQ:PL  1/2:0,32,12:55:99:4008,964,1808,3176,0,5166")
        assert v.equalsAnnovar("chr3", 40503526, "CTGCTGCTG")
        assert v.equalsAnnovar("chr3", 40503521, "-")
        
        v = Variant.parse("chr3 46751073    rs10578999  TAAGAAG T,TAAG  32852.19    PASS    AC=1,1;AF=0.500,0.500;AN=2;DB;DP=158;FS=0.000;HaplotypeScore=84.0101;MLEAC=1,1;MLEAF=0.500,0.500;MQ=59.01;MQ0=0;QD=207.93;RPA=9,7,8;RU=AAG;STR;set=Intersection GT:AD:DP:GQ:PL  1/2:0,12,125:155:99:6991,6090,6305,609,0,197")
        assert v.equalsAnnovar("chr3", 46751077, "-")
        
        v = Variant.parse("chr3 46751077    rs10578999  A T,C  32852.19    PASS    AC=1,1;AF=0.500,0.500;AN=2;DB;DP=158;FS=0.000;HaplotypeScore=84.0101;MLEAC=1,1;MLEAF=0.500,0.500;MQ=59.01;MQ0=0;QD=207.93;RPA=9,7,8;RU=AAG;STR;set=Intersection GT:AD:DP:GQ:PL  1/2:0,12,125:155:99:6991,6090,6305,609,0,197")
        assert v.equalsAnnovar("chr3", 46751077, "C") == 2
        
        v = Variant.parse("chr3 46751077    rs10578999  A T,C  32852.19    PASS    AC=1,1;AF=0.500,0.500;AN=2;DB;DP=158;FS=0.000;HaplotypeScore=84.0101;MLEAC=1,1;MLEAF=0.500,0.500;MQ=59.01;MQ0=0;QD=207.93;RPA=9,7,8;RU=AAG;STR;set=Intersection GT:AD:DP:GQ:PL  1/2:0,12,125:155:99:6991,6090,6305,609,0,197")
        assert v.equalsAnnovar("chr3", 46751077, "T") == 1
        
        v = Variant.parse("chr4 147560457   rs5862765   TGGCGGCGGCGGC   TGGCGGCGGCGGCGGC,TGGC,T 480.14  .   AC=29,6,2;AF=0.604,0.125,0.042;AN=48;DB;DP=470;MQ0=0;RU=GGC;STR;set=variant-variant2-variant10-variant11-variant13-variant15-variant17-variant19-variant20-variant21-variant22-variant24-variant25-variant26-variant28-variant30-variant31-variant42-variant43-variant45-variant46-variant52-variant54-variant57    GT:DP:GQ    0/1:21:10   0/2:14:61   ./. ./. ./. ./. ./. ./. ./. 1/1:9:15    1/1:13:18   ./. 0/2:24:99   ./. 1/1:12:30   ./. 1/1:36:25   ./. 0/1:17:99   1/1:19:21   0/2:18:99   0/1:27:93   ./. 1/3:27:99   0/3:20:99   1/1:17:21   ./. 1/1:18:15   ./. 1/1:17:21   1/1:22:24   ./. ./. ./. ./. ./. ./. ./. ./. ./. ./. 1/1:10:18   1/1:22:24   ./. 0/1:26:68   0/2:17:99   ./. ./. ./. ./. ./. 0/2:16:99   ./. 0/2:17:99   ./. ./. 1/1:14:24   ./.")
        assert v.equalsAnnovar("chr4", 147560469, "GGC") == 1
    }

}
