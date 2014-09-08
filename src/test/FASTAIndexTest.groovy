import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.Test;


class FASTAIndexTest {
    
    def seq1 = "GTCACCAGTACCAGGGACAGAGCCCAGGTGGGGTGGGGGCGGGGTCCAGCACCACGGCCAGCACCGACCACCAGGACCCCGGAGCCAGCACCATGGACAGAAAACTGCCCACCAGGATCTGACGCCAGCACGCCGCCAGGCCCACACAGGGTCTCCGGTCAGAGTCCCAGGGTCAGCTCCCAGCAGGGCCTAGGGGAGGCTGGACCAGCTCCCTGTGCCTCATTCCAAGGCAGCCCAGCCGGAGAGAAGGGGCACAGGCCACACATCTGTCCCATAAAATTAAACGCTTTTTAGTGTTTAAAATAAGCAGCATTTACACAGAAGCAGCTCTATGTTAACCATCTAAACGCTGGGACTTTGATACA"
    def seq1_name = "chr1:1270388-1270753"
    
    def seq2 = "GTCAGCTCCCAGCAGGGCCTAGGGGAGGCTGGACCAGCTCCCTGTGCCTCATTCCAAGGCAGCCCAGCCGGAGAGAAGGGGCACAGGCCACACATCTGTCCCATAAAATTAAACGCTTTTTAGTGTTTAAAATAAGCAGCATTTACACAGAAGCAGCTCTATGTTAACCATCTAAA"    
    def seq2_name = "chr1:1270559-1270735"
    
    FASTAIndex index 
    
    @Before
    void setup() {
        index = new FASTAIndex(new FASTA("tests/test.fasta"))
    }
    
    @Test
    public void testIndex() {
        
        println index.querySequenceName(seq1) 
        
//        println index.sequences
        
        // Basic lookup of a sequence should work
        assert index.querySequenceName(seq1) == seq1_name
        assert index.querySequenceName(seq2) == seq2_name
        
        // Lookup with missing 2 chars of prefix should work
        assert index.querySequenceName(seq1.substring(2)) == seq1_name
        assert index.querySequenceName(seq2.substring(2)) == seq2_name
        
        // Won't find it at large offset
        assert index.querySequenceName(seq1.substring(10)) == null
    }
    
    @Test
    public void testLargeOffset() {
        
        index = new FASTAIndex(new FASTA("tests/test.fasta"), 0..20)
        
        println "Querying for subsequence ${seq1.substring(10)} from sequences: ${index.sequences.keySet().join('\n')}"
        
        println index.sequences[seq1.substring(10).substring(0,30)]
        
        assert index.querySequenceName(seq1.substring(10)) == seq1_name
    }
    
    @Test
    public void testSmallMax() {
        index = new FASTAIndex(new FASTA("tests/test.fasta"), 0..5, 2)
        
        // Should not index the second sequence because it is over the max limit
        assert index.querySequenceName(seq1) == seq1_name
        assert index.querySequenceName(seq2) == null
    }
    
    @Test
    public void testFindAmplicon() {
        
        println "Loading regions ..."
        BED regions = new BED("/Users/simon/work/dsd/batch7.bowtie/amplicons.chr1.bed").load()
        
        println "Indexing sequences ..."
        FASTAIndex index = new FASTAIndex(new FASTA("/Users/simon//hg19/ucsc.hg19.fasta"), regions)
        
        // Try to find the amplicon for a known sequence
        def result = index.querySequenceName("CTTGGAATCCCGGTCATCAGACTAAGAAGCAAAGAGAAAGTTAATCATTTTCTTTATAAGTTTTTTTTTCTTCATAGATAAAAGTATTTTTAATGATAATCAAACCTGGGCAACATCCCAAAACAAACTTTCACATGTACTCTGAATGAGC")
        
        println "Found sequence $result"
        
    }
}
