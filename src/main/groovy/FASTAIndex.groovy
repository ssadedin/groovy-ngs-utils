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

import groovy.transform.CompileStatic;
import net.sf.samtools.util.SequenceUtil;

/**
 * A very simplistic FASTA index that allows lookup of sequence names by their
 * sequence content using a fixed length prefix seed. It only supports looking up
 * sequences by their prefixes, with an optional offset
 * from the query start by up to specified number of bases (default=5). Both 
 * the original sequence AND its reverse complement are indexed, so you can
 * perform a single query to identify a sequence in the case you are not sure
 * what the strand / orientation of the query sequence is.
 * <p>
 * NOTE: this class is intended for indexing large numbers of SHORT sequences.
 * It will not work for indexing, for example, a reference sequence for an 
 * organism! (You will be able to look up each chromosome by a short prefix, not
 * terribly useful).
 * <p>
 * NOTE2: this class is not for working at low level with pre-indexed FASTA files 
 * (eg: .fai format) For working with those, use the {@link FASTA} class. This class
 * takes a {@link FASTA} object and adds ability to look up by sequence content
 * to it.
 * <p>Example:
 * <pre>
 * index = new FASTAIndex(new FASTA("tests/test.fasta"), 0..20)
 * assert index.querySequence("AGTCCCTATTACAAA") == "amplicon_1"
 * </pre>
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class FASTAIndex {

    /**
     * Index of amplicon names maps name => full sequence
     */
    Map<String,String> sequenceNames = new HashMap(40000)
    
    /**
     * Index of amplicon sequences, maps subsequence to amplicon name(s)
     */
    Map<String,String> sequences = new HashMap(40000)
    
    /**
     * Range of offsets from beginning of sequences to index
     * (memory expensive to increase this a lot)
     */
    IntRange offsetRange = 0..6
    
    /**
     * Maximum number of sequences to index (0 means unlimited)
     */
    int maxSize = 0 
    
    /**
     * Size of seed to use (bp)
     */
    int seedSize = 30
    
    /**
     * Create an index from the given fasta, where each fasta sequence corresponds
     * to a single amplicon
     * 
     * @param fasta
     */
    public FASTAIndex(FASTA fasta, Regions regions=null) {    
        this.index(fasta,regions)
    }
	
	/**
	 * For unit tests only
	 */
	public FASTAIndex() {
	}
    
    /**
     * Index the given fasta, masked using the given BED file
     * @param fasta
     * @param bed
     */
    void index(FASTA fasta, Regions regions) {
        ProgressCounter counter = new ProgressCounter()
        
        int countTooShort = 0
        
        def handler = { String amplicon, byte[] bases ->
            
            // Ignore stupidly short sequences
            if(bases.size()<seedSize+offsetRange.to+1) {
                ++countTooShort
                return
            }
            
            boolean log = false
//            if(amplicon=="chr1:154926092-154926234") {
//                log=true
//            }
            
            sequenceNames[amplicon] = new String(bases).toUpperCase()
            
            for(i in offsetRange.from..(offsetRange.to+1)) {
              String sequence = new String(bases,i,seedSize).toUpperCase()
              if(log)
                  println "Indexing $amplicon: $sequence"
              sequences[sequence] = amplicon
            }
            
            SequenceUtil.reverseComplement(bases)
            for(i in offsetRange.from..(offsetRange.to+1)) {
              String sequence = new String(bases,i,seedSize).toUpperCase()
              sequences[sequence] = amplicon
            }
            
            counter.count()
            
            if(maxSize > 0 && counter.count >= maxSize) 
                throw new Abort()
        }
        
        if(regions) {
           fasta.eachSequence(regions, handler)
        }
        else {
           fasta.eachSequence(handler)
        }
        
        if(countTooShort * 100 / counter.count > 10) {
            throw new IllegalStateException("More than 10% of your sequences were shorter than the seed length. Please use a shorter seed length or check that your data is correct.")
        }
    }
    
    String querySequenceName(String sequence) {
        if(sequence.size()<seedSize) 
            System.err.println("A sequence queried ($sequence) was shorter than the seed length ($seedSize). Querying sequences shorter than the seed size is not supported.")
        return sequences[sequence.substring(0, seedSize)]
    }
    
    /**
     * Create an index from the given fasta, where each fasta sequence corresponds
     * to a single amplicon
     * 
     * @param fasta
     */
    public FASTAIndex(FASTA fasta, IntRange offsetRange, int maxSize=0, int seedSize=30, BED bed=null) {
        this.maxSize = maxSize
        this.seedSize = seedSize
        this.offsetRange = offsetRange
        this.index(fasta,bed)
    }
}