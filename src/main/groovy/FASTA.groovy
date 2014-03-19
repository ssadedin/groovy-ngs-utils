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
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;


/**
 * Simple wrapper to support groovy interface & utilities for FASTA files
 * <p>
 * Only indexed files are supported
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class FASTA {
    
    public static final long PRINT_INTERVAL_MS = 15000
    
    IndexedFastaSequenceFile indexedFastaFile
    
    FASTA(String fastaFile) {
        this(new File(fastaFile))
    }
    
    FASTA(File fastaFile) {
        this.indexedFastaFile  = new IndexedFastaSequenceFile(fastaFile)
    }
    
    /**
     * Returns the sequence of bases over the given range
     * <b>NOTE:</b>The range is <i>inclusive</i>.
     * 
     * @param contig    contig to query
     * @param start     start position (inclusive)
     * @param end       end position (inclusive)
     * @return
     */
    String basesAt(String contig, long start, long end) {
      return new String(this.indexedFastaFile.getSubsequenceAt(contig, start, end).bases)
    }
    
    byte[] baseBytesAt(String contig, long start, long end) {
      return this.indexedFastaFile.getSubsequenceAt(contig, start, end).bases
    }
     
    public static final int T = (int)"T".charAt(0)
    public static final int A = (int)"A".charAt(0)
    public static final int C = (int)"C".charAt(0)
    public static final int G = (int)"G".charAt(0)
    public static final int R = (int)"R".charAt(0)
    public static final int Y = (int)"Y".charAt(0)
    public static final int W = (int)"W".charAt(0)
    public static final int K = (int)"K".charAt(0)
    public static final int M = (int)"M".charAt(0)
    public static final int S = (int)"S".charAt(0)
    public static final int H = (int)"H".charAt(0)
    public static final int D = (int)"D".charAt(0)
    public static final int B = (int)"B".charAt(0)
    public static final int N = (int)"N".charAt(0)
    public static final int V = (int)"V".charAt(0)
    
    @CompileStatic
    static String reverseComplement(String bases) {
        byte [] bytes = bases.bytes
        StringBuilder result  = new StringBuilder()
        for(int i=bytes.length-1; i>=0; --i) {
            switch(bytes[i]) {
                case A:
                    result.append((char)T)
                    break
                case T:
                    result.append((char)A)
                    break
                case C:
                    result.append((char)G)
                    break
                case G:
                    result.append((char)C)
                    break
                case Y:
                    result.append((char)R)
                    break
                case R:
                    result.append((char)Y)
                    break
                case W:
                    result.append((char)W)
                    break
                case M:
                    result.append((char)K)
                    break
                case K:
                    result.append((char)M)
                    break
                case S:
                    result.append((char)S)
                    break
                case H:
                    result.append((char)D)
                    break
                case D:
                    result.append((char)H)
                    break
                case B:
                    result.append((char)V)
                    break
                case V:
                    result.append((char)B)
                    break
                case N:
                    result.append((char)N)
                    break
            }
        }
        return result.toString()
    }
    
    /**
     * Iterates over the FASTA and calls the given closure for each
     * reference sequence (eg: chromosome) passing the bases of the sequence
     * as a byte array
     */
    @CompileStatic
    void eachSequence(Closure c) {
        ReferenceSequence seq
        int count=0
        long lastPrintTimeMs = System.currentTimeMillis()
        while((seq = this.indexedFastaFile.nextSequence()) != null) {
            c(seq.name, seq.bases)
            ++count
            if((count % 10000 == 0) && (System.currentTimeMillis() - lastPrintTimeMs) > PRINT_INTERVAL_MS) {
                System.err.println("Processed " + count + " reads")
                lastPrintTimeMs = System.currentTimeMillis()
            }
        }
    }
    
    /**
     * Iterates over the ranges in teh BED file and extracts the bases for the ranges,
     * passing the bases of the sequence
     * as a byte array to the given closure for each range.
     */
    @CompileStatic
    void eachSequence(BED bed, Closure c) {
        boolean includeRegion = false
//        if(c.maximumNumberOfParameters > 2) 
//            includeRegion = true
        bed.eachRange(unique:true) {  String chr, int start, int end ->
            c(chr+":"+start+"-"+end, this.baseBytesAt(chr,start,end))
        }
    }
    
    @CompileStatic
    static String format(String contig, String sequence) {
        StringBuilder output = new StringBuilder()
        final int iterations = (int)Math.ceil(sequence.size()/80)
        for(int i=0; i<iterations; ++i) {
            output.append(sequence.substring(i*80, Math.min((i+1)*80,sequence.size())))
            output.append("\n")
        }
        return output.toString()
    }
    
    /**
     * Display given sequence with a single base bracketed to highlight it
     * 
     * @param seq   sequence to display
     * @param start position to start from (inclusive)
     * @return  sequence with specified base bracketed
     */    
    static String bracket(String seq, int position) {
        bracket(seq,position,position+1)
    }
    
    /**
     * Display given sequence with given bases bracketed
     * 
     * @param seq   sequence to display
     * @param start position to start from (inclusive)
     * @param end   position to end at (exclusive)
     * @return  sequence with specified bases bracketed
     */
    static String bracket(String seq, int start, int end) {
        "${seq.substring(0,start)}[${seq.substring(start,end)}]${seq.substring(end)}"
    }
}
