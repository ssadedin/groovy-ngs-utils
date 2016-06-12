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

import htsjdk.samtools.util.SequenceUtil;

import org.junit.Test;

import RestrictionEnzyme as RE

class RestrictionEnzymeTest {
    
    RE re

    @Test
    public void testNoWildCards() {
        
        re = new RE("x", "CATG", 2,0)
        
        assert re.cutsSequenceBefore("CATG", 2)
        assert !re.cutsSequenceBefore("CATG", 1)
        assert !re.cutsSequenceBefore("CATG", 3)
        assert !re.cutsSequenceBefore("CATG", 4)
        assert !re.cutsSequenceBefore("CATG", 0)
    }

    @Test
    public void testSingleWildCard() {
        
        def re = new RE("x", "CANTG", 2,0)
        
        assert re.cutsSequenceBefore("CAATG", 2)
        assert re.cutsSequenceBefore("CACTG", 2)
        assert re.cutsSequenceBefore("CATTG", 2)
        assert re.cutsSequenceBefore("CAGTG", 2)
        
        assert !re.cutsSequenceBefore("CATG", 1)
        assert !re.cutsSequenceBefore("CTATG", 2)
        assert !re.cutsSequenceBefore("CTATG", 3)
        assert !re.cutsSequenceBefore("CAAGG", 2)
        assert !re.cutsSequenceBefore("CAAGG", 3)
        assert !re.cutsSequenceBefore("CAAGG", 4)
    }
    
    @Test
    public void testTrailingWildCards() {
        
        re = new RE("x", "CANNNNN", 2,2)
        
        assert re.cutsSequenceBefore("CAATGAA", 2)
        assert re.cutsSequenceBefore("CATGGGG", 2)
        
        assert !re.cutsSequenceBefore("GATGGGG", 2)
        assert !re.cutsSequenceBefore("CGTGGGG", 2)
    }
    
    @Test
    public void testLongSequence() {
        re = new RE("x", "CCATCNNNNN",8, 1)
        
        // This string has GGCC at 18,19,20,21
        def seq = "GGGTAAGGGTTACAGGTGGGCCTTAGAAGGAGCCAGGGGG" 
        
        assert !re.cutsSequenceBefore(seq, 20)
        
        re = new RE("x", "GGCC",2, 0)
        assert re.cutsSequenceBefore(seq, 20)
        
        re = new RE("x", "AATT",0, 4)
        assert !re.cutsSequenceBefore(seq, 20)
        
        assert re.cutsSequenceBefore("AATT", 4)
        assert re.cutsSequenceBefore("AATT", 0)
        
        re = new RE("x", "CATG", 4, -4)
        assert !re.cutsSequenceBefore(seq, 20)
    }
    
    
    @Test
    public void testOverHangCut() {
        re = new RE("x", "CCTT",2, 0)
        
        assert !re.cutsSequenceBefore("TTCC", 2)
        assert re.cutsSequenceBefore("CCTT", 2)
        
        re = new RE("x", "CCTT",1, 2)
        
        assert re.cutsSequenceBefore("CCTT", 3)
        
        re = new RE("x", "CCATCNNNNNN", 9, 1)
        assert !re.cutsSequenceBefore("TTCCCCCTGCATCTTCCTCAGACTTGCAGATGATTTTACT", 20)
    }
    
    @Test
    public void testReverseComplement() {
        
        def re = new RE("x", "CCTCNNNNNNN", 11, -1)
        println "Reverse complement = " + SequenceUtil.reverseComplement("CCTCNNNNNNN")
        assert re.cutsSequenceBefore("CCCAGGAGAGTGTAGAAGGGGCAGGGGGAGGGAGCGTGTC", 20)
    }
    
    @Test
    public void testDebugMspJI() {
        def re = new RE("x", "CNNRNNNNNNNNNNNNN",13,4);
        
        assert re.cutsSequenceBefore("GGGGTGAACGGAGTCGCCGCTTTTTAAACAGCCTGGGGTC", 20)
    }
    
    @Test
    public void testDebugBsp1286I() {
        re = new RE("x", "GDGCHC",5,-4)
        assert re.cutsSequenceBefore("GCAGCCCTGGTCGCAGTGCCCGTCGCTGAAGTGGTCCTTG", 20)
    }
    
    @Test
    void testIUPAC() {
        def seq = "ACAATGGAGCTATCTCATCTAAGTGATTCCTGATGCCAAA"
        RE re = new RE("x", "CAYNNNNRTG",5,0)
        assert re.cutsSequenceBefore(seq, 20)
        
        // Now what about reverse complement?
        println "Reverse complement of CAYNNNNRTG = " + SequenceUtil.reverseComplement("CAYNNNNRTG")
        println "Reverse complement of CAYNNNNRTG = " + FASTA.reverseComplement("CAYNNNNRTG")
        assert re.cutsSequenceBefore(SequenceUtil.reverseComplement(seq), 20)
    }
}
