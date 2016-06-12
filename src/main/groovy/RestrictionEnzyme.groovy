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

import java.util.regex.Pattern;

import htsjdk.samtools.util.SequenceUtil;
import static FASTA.*

/**
 * A simple class for matching motifs in restriction enzymes to 
 * FASTA sequence. Only works with cutters that recognise 
 * fixed length motifs for now.
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class RestrictionEnzyme {
    
    String name
    Pattern forwardSuffixPattern
    Pattern forwardPrefixPattern
    Pattern reversePrefixPattern
    Pattern reverseSuffixPattern
    
    RestrictionEnzyme complement
    
    String motif
    int motifSize
//    List<Range> recognitionSites
    
    String forwardSuffix
    String forwardPrefix
    String reversePrefix
    String reverseSuffix
   
    /**
     * Create a restriction enzyme object with given name, motif cut point and overhang.
     */
    RestrictionEnzyme(String name, String motif, int forwardCutPoint, int overHang, boolean withComplement=true) {
        
        this.name = name
        this.motifSize = motif.size()
        this.motif = motif
        
        if(forwardCutPoint > motifSize || forwardCutPoint+ overHang>motifSize)
            throw new IllegalArgumentException("Enzyme cut points beyond the length of specified motif")
        
        // Motif after cut point
        if(forwardCutPoint<motifSize) {
            forwardSuffix = motif.substring(forwardCutPoint)
            def p = toRegex(forwardSuffix)
            this.forwardSuffixPattern = Pattern.compile(p)
        }
        
        // Motif before cut point
        if(forwardCutPoint>0) {
            forwardPrefix = motif.substring(0,forwardCutPoint)
            def p = toRegex(forwardPrefix)
            this.forwardPrefixPattern = Pattern.compile(p)
        }
        
        int reverseCutPoint = forwardCutPoint + overHang
        if(reverseCutPoint > 0 && reverseCutPoint<=motifSize) {
            reversePrefix = motif.substring(0, reverseCutPoint)
            def p = toRegex(reversePrefix)
            this.reversePrefixPattern = Pattern.compile(p)
        }
        
        // Motif before cut point
        if(reverseCutPoint>=0) {
            reverseSuffix = motif.substring(reverseCutPoint)
            def p = toRegex(reverseSuffix)
            this.reverseSuffixPattern = Pattern.compile(p)
        }
        
        if(withComplement) 
            this.complement = new RestrictionEnzyme(name, FASTA.reverseComplement(motif), motifSize - forwardCutPoint, -1*overHang, false)
    }
    
    /**
     * Convert N's to wildcards and select IUPAC codes to appropriate character sets
     */
    String toRegex(String motif) {
        
        char [] motifChars = motif.toCharArray()
        
        int NOT_REC = 0, REC=1;
        int state = NOT_REC;
        Range currentRec = null;
        StringBuilder result = new StringBuilder()
        
//        this.recognitionSites = []
        
        for(int i=0; i<motifChars.length; ++i) {
            switch((int)motifChars[i]) {
                case N:
//                    if(state == REC) {
//                        state = NOT_REC
//                        currentRec.to = i
//                        this.recognitionSites << currentRec
//                        currentRec = null
//                    }
                    result.append(".")
                    break
                    
                case Y:
                    result.append("[CT]")
                    break
                case R:
                    result.append("[AG]")
                    break
                case W:
                    result.append("[AT]")
                case D:
                    result.append("[AGT]")
                    break
                case H:
                    result.append("[ACT]")
                    break
                default:
                    result.append((char)motifChars[i])
            }
//            
//            if(motifChars[i] != N && state == NOT_REC) {
//                currentRec = new IntRange(i,-1)
//                state = REC
//            }
            
        }
        return result.toString()
    }
    
    /**
     * Returns true if this restriction enzyme would cut the given sequence
     * in between the given position and the previous position.
     * 
     * @param sequence 
     * @param position zero-based position. ie: 0 means the cut would be prior to the first base in the sequence
     * @return  true iff the enzyme would cut the specified position
     */
    boolean cutsSequenceBefore(String sequence, int position) {
        
        if(checkMatchForward(sequence,position) ||  checkMatchReverse(sequence,position))
            return true
            
        if(this.complement != null)
            return this.complement.cutsSequenceBefore(sequence, position)
    }
    
    
    private boolean checkMatchReverse(String fasta, int position) {
        if(reversePrefixPattern != null) {
            
            if(position<reversePrefix.size())
                return false
            
            String reversePrefix = fasta.substring(position - reversePrefix.size(), position)
            if(!reversePrefixPattern.matcher(reversePrefix).matches())
                 return false
//            println "Match fasta $fasta at position $position with motif $reversePrefixMotif.pattern prefix = $reversePrefix"
        }
        
        if(reverseSuffixPattern != null) {
            int endIndex = position + reverseSuffix.size()
            if(endIndex > fasta.size())
                return false
                
            String reverseSuffix = fasta.substring(position, endIndex)
            if(!reverseSuffixPattern.matcher(reverseSuffix).matches()) {
                 return false
            }
        }
         
        // Both match so actual pattern matches too
        return true    
        
    }
    
    private boolean checkMatchForward(String fasta, int position) {
        
        if(forwardPrefixPattern != null) {
            
            int patternSize = forwardPrefix.size()
            if(position < patternSize)
                return false
                
            String prefix = fasta.subSequence(position - patternSize, position)
            if(!forwardPrefixPattern.matcher(prefix).matches())
                return false
        }
        
        
        if(forwardSuffixPattern != null) {
            int endIndex = position+forwardSuffix.size()
            if(endIndex>fasta.size())
                return false
            String suffix = fasta.substring(position,endIndex)
            if(!forwardSuffixPattern.matcher(suffix).matches())
                return false
        }
        
        return true
    }
    
    String toString() {
        "$name[${forwardPrefixPattern?.pattern?:''}^${forwardSuffixPattern?.pattern?:''}]"
    }
    
    public static List load(String fileName) {
        load(new FileInputStream(fileName))
    }
    
    /**
     * Convenience method to load enzymes from UCSC's cutter table format
     * 
     * @param inStream  Input stream to load from
     * @return  list of enzymes loaded
     */
    public static List load(InputStream inStream) {
        // For each enzyme, we generate 2 patterns
        // The forward pattern must match AFTER the cut point
        // The reverse pattern must match leading UP TO the cut point
        List<RestrictionEnzyme> enzymes = []
        
        // Read our file of cutters
        inStream.eachLine { line ->
            if(line.startsWith("#"))
                return
            def fields = line.split(/\t/)
            def enzyme = new RestrictionEnzyme(fields[0], fields[3], fields[4].toInteger(), fields[5].toInteger())
            enzymes << enzyme
        }
        
        return enzymes
    }
    
    /**
     * Enzymes possibly used by HaloPlex, extracted from the UCSC cutters database
     */
    final static String HALOPLEX_DEFAULT_ENZYMES = """
        BccI    10  5   CCATCNNNNN  9   1   0   0   0       1   N,  3   571,576,585,
        DdeI    5   4   CTNAG   1   3   1   0   2   BstDEI,HpyF3I,  7   M,N,O,Q,R,S,X,  4   79,263,316,517,
        HaeIII  4   4   GGCC    2   0   1   0   4   BshFI,BsnI,BsuRI,PhoI,  13  B,I,J,K,M,N,O,Q,R,S,U,X,Y,  3   77,527,544,
        HpyCH4V 4   4   TGCA    2   0   1   0   0       1   N,  1   588,
        MnlI    11  4   CCTCNNNNNNN 11  -1  0   0   0       6   F,I,N,Q,V,X,    5   76,303,435,748,998,
        AluI    4   4   AGCT    2   0   1   0   1   AluBI,  16  B,C,F,I,J,K,M,N,O,Q,R,S,U,V,X,Y,    6   432,447,727,796,948,993,
        Tsp509I 4   4   AATT    0   4   1   0   4   TspEI,MluCI,Sse9I,TasI, 1   N,  2   585,845,
        MseI    4   4   TTAA    1   2   1   0   3   SaqAI,Tru1I,Tru9I,  2   B,N,    3   577,866,867,
        HpyCH4III   5   4   ACNGT   3   -1  1   0   2   Bst4CI,TaaI,    1   N,  2   590,1020,
        MspJI   17  2   CNNRNNNNNNNNNNNNN   13  4   0   0   0       1   N,  3   162,1009,1013,
        NlaIII  4   4   CATG    4   -4  1   0   5   CviAII,FaeI,FatI,Hin1II,Hsp92II,    1   N,  3   446,578,690,
        BsaJI   6   4   CCNNGG  1   4   1   0   2   BseDI,BssECI,   1   N,  3   416,585,977,
    """
}