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

import net.sf.samtools.util.SequenceUtil;
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
}