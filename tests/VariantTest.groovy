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
    
    @Test void testFixGenotypeOrder() {
        v = new Variant(ref:"CAG", alt:"C", pos: 20)
        v.line = "chr1  57179502    .   G   T   51.8    PASS    AC=1;AF=0.5;AN=2;BaseQRankSum=-4.147;DP=42;Dels=0;FS=17.681;HaplotypeScore=72.6709;MLEAC=1;MLEAF=0.5;MQ=41.17;MQ0=0;MQRankSum=-2.887;QD=1.23;ReadPosRankSum=4.331;set=Intersection;EFF=DOWNSTREAM(MODIFIER||||728|C1orf168||CODING|NM_001004303|),UTR_3_PRIME(MODIFIER||||552|PRKAA2||CODING|NM_006252|9)   AD:DP:GQ:GT:PL  .:40:80:0/1:80,0,255"
        v.fixGenoTypeOrder()
        
        println v.line
        
        assert v.line.split('\t')[8] == 'GT:AD:DP:GQ:PL'
        assert v.line.split('\t')[9] == '0/1:.:40:80:80,0,255'
    }

}
