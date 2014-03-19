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
 * sequences. Currently indexes the first 30 base pairs of each sequence,
 * and that is all!
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class FASTAIndex {

    /**
     * Index of amplicon names maps name => sequence
     */
    Map sequenceNames = new HashMap(40000)
    
    /**
     * Index of amplicon sequences, maps subsequence to amplicon name(s)
     */
    Map sequences = new HashMap(40000)
    
    /**
     * Create an index from the given fasta, where each fasta sequence corresponds
     * to a single amplicon
     * 
     * @param fasta
     */
    public FASTAIndex(FASTA fasta, BED bed=null) {
        ProgressCounter counter = new ProgressCounter()
        
        def handler = { String amplicon, byte[] bases ->
            
            // Ignore stupidly short sequences
            if(bases.size()<40)
                return
            
            sequenceNames[amplicon] = new String(bases).toUpperCase()
            
            for(i in 0..6) {
              String sequence = new String(bases,i,30).toUpperCase()
              sequences[sequence] = amplicon
            }
            
            SequenceUtil.reverseComplement(bases)
            for(i in 0..6) {
              String sequence = new String(bases,i,30).toUpperCase()
              sequences[sequence] = amplicon
            }
            counter.count()
        }
        
        if(bed) {
           fasta.eachSequence(bed, handler)
        }
        else {
           fasta.eachSequence(handler)
        }
    }
}