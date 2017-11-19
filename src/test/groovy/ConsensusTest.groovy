/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2016 Simon Sadedin, ssadedin<at>gmail.com
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

import gngs.Utils

class ConsensusTest {
    
    @Test
    void testOneSequence() {
        
        def cons = new Consensus(["ATCGCGA"])
        cons.build()
        assert cons.bases == "ATCGCGA"
        assert cons.score == 0.0d
    }
    
    @Test
    void testOneDifference() {
        def cons = new Consensus([
            "ATCGCGA",
            "ATCGCGA",
            "ATCGCGT",
        ])
        cons.build()
        assert cons.bases == "ATCGCGA"
        assert cons.score > 6.0d  && cons.score < 7.0d
    }

    @Test
    void testTwoDifferences() {
        def cons = new Consensus([
            "ATCGCGA",
            "ATCGCGA",
            "ATTGCGA",
            "ATCGCGT",
        ])
        cons.build()
        assert cons.bases == "ATCGCGA"
        assert cons.score > 6.0d  && cons.score < 7.0d
    }
    
    @Test
    void testTwoDifferencesSameBase() {
        def cons = new Consensus([
            "ATCGCGA",
            "ATCGCGA",
            "ATCGCGT",
            "ATCGCGT",
        ])
        cons.build()
        assert cons.bases == "ATCGCGA"
        assert cons.score > 6.0d  && cons.score < 6.4d
    }
    
    
    @Test
    void testOppositeOrder() {
        def cons = new Consensus([
            "ATCGCGA",
            "ATCGCGA",
            "ATTGCGA",
            "ATCGCGT",
        ].reverse())
        cons.build()
        assert cons.bases == "ATCGCGA"
        assert cons.score > 6.0d  && cons.score < 7.0d
    } 
    
    @Test
    void testAllDifferent() {
        def cons = new Consensus([
            "ATCGCGA",
            "TCACGCT",
            "CGTTATG",
            "GAGATAC",
        ].reverse())
        cons.build()
        assert cons.bases == "AAAAAAA"
        assert cons.score == 0.0d
    }  
    
    @Test
    void testRead() {
        def cons = new Consensus([
            "CTGTCCCTGCCCCACCTCCACCCCCGCCACCTCCCCCACCTGCCACCCCTGTGACCCCGGCCCCCGTGCCTCCCTTCGAGAAGCAAGGAGGAAAGGACAAGGAAGACAAGCAGACATTCCAAGTCACAGACTGTCGAAGTTTGGTCAAAAC",
            "CTGTCCCTGCCCCACCTCCACCCCCGCCCCCACCCCCACCTGCCACCCCTGTGACCCCGGCCCCCGTGCCTCCCTTCGAGAAGCAAGGAGAAAAGGACAAGGAAGACAAGCAGACATTCCAAGTCACAGACTGTCGAAGTTTGGTCAAAAC",
            "CTGTCCCTGCCCCACCTCCCCCCCCGCCCCCACCCCCACCTGCCACCCCTGTGACCCCGGCCCCCGTGCCTCCCTTCGAGAAGCAAGGAGAAAAGGACAAGGAAGACAAGCAGACATTCCAAGTCACAGACTGTCGAAGTTTGG",
        ])
        cons.build()
        assert cons.bases.startsWith("CTGTCCCTGCCCCACCTCCACCCCCGCCCCCACCCCCACCTGCCACCCCTGTGACCCCGGCCCCCGTGCCTCCCTTCGAGAAGCAAGGAGAAAAGGAC")
        assert cons.score > 10.0d
    }
    
    @Test
    void testPerformance() {
        
        long timeMs = Utils.timeMs {
            for(i in 0..1000) {
                def cons = new Consensus([
                    "ATCGCGA"*2,
                    "ATCGCGA"*2,
                    "ATCGCGT"*2,
                    "ATCGCGT"*2,
                ])
                cons.build()
            }
        }
        
        assert timeMs < 10000
        
        println "Average time = " + (timeMs / 1000.0d)
    }
}
