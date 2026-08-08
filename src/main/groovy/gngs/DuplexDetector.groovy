package gngs

import groovy.transform.CompileStatic
import htsjdk.samtools.Cigar
import htsjdk.samtools.CigarElement
import htsjdk.samtools.CigarOperator
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.util.SequenceUtil as SAMSequenceUtil

/**
 * Detects Oxford Nanopore reads that were accidentally read in duplex form
 * but processed as simplex reads.
 * 
 * A duplex read has:
 * - An aligned part (mapped to reference)
 * - A soft clipped part (unmapped, at the end of the read)
 * - Both parts are approximately the same length (within tolerance)
 * - The soft clipped part is the reverse complement of the aligned part
 * 
 * @author simon.sadedin@mcri.edu.au
 */
@CompileStatic
class DuplexDetector {
    
    /**
     * Default tolerance for length comparison (10%)
     */
    static final double DEFAULT_LENGTH_TOLERANCE = 0.2
    
    /**
     * Default alignment threshold
     */
    static final double DEFAULT_ALIGNMENT_THRESHOLD = 0.30
    
    /**
     * Default number of bases to compare for alignment check
     */
    static final int DEFAULT_ALIGNMENT_WINDOW_SIZE = 50
    
    /**
     * Number of bases to compare for alignment check
     */
    int alignmentWindowSize = DEFAULT_ALIGNMENT_WINDOW_SIZE
    
    /**
     * Tolerance for comparing aligned vs soft clipped lengths
     */
    double lengthTolerance = DEFAULT_LENGTH_TOLERANCE
    
    /**
     * Threshold for considering sequences aligned (fraction of matching bases)
     */
    double alignmentThreshold = DEFAULT_ALIGNMENT_THRESHOLD
    
    /**
     * Information about soft clips in a read
     */
    static class SoftClipInfo {
        int leadingSoftClip = 0
        int trailingSoftClip = 0
        int alignedBases = 0
        
        int getTotalSoftClip() {
            return leadingSoftClip + trailingSoftClip
        }
        
        boolean hasSoftClips() {
            return leadingSoftClip > 0 || trailingSoftClip > 0
        }
    }
    
    /**
     * Extract soft clip information from a CIGAR string
     * 
     * @param record SAMRecord to analyze
     * @return SoftClipInfo containing counts of leading, trailing, and aligned bases
     */
    SoftClipInfo extractSoftClipInfo(SAMRecord record) {
        SoftClipInfo info = new SoftClipInfo()
        
        Cigar cigar = record.getCigar()
        if (cigar == null || cigar.numCigarElements() == 0) {
            return info
        }
        
        List<CigarElement> elements = cigar.getCigarElements()
        
        // Check for leading soft clip (skip any leading hard clips)
        for (CigarElement element : elements) {
            if (element.getOperator() == CigarOperator.SOFT_CLIP) {
                info.leadingSoftClip = element.getLength()
                break
            } else if (element.getOperator() != CigarOperator.HARD_CLIP) {
                // Stop if we hit something other than H or S
                break
            }
        }
        
        // Check for trailing soft clip (skip any trailing hard clips)
        for (int i = elements.size() - 1; i >= 0; i--) {
            CigarElement element = elements[i]
            if (element.getOperator() == CigarOperator.SOFT_CLIP) {
                info.trailingSoftClip = element.getLength()
                break
            } else if (element.getOperator() != CigarOperator.HARD_CLIP) {
                // Stop if we hit something other than H or S
                break
            }
        }
        
        // Count aligned bases (M, =, X operators)
        for (CigarElement element : elements) {
            CigarOperator op = element.getOperator()
            if (op == CigarOperator.M || op == CigarOperator.EQ || op == CigarOperator.X) {
                info.alignedBases += element.getLength()
            }
        }
        
        return info
    }
    
    /**
     * Check if soft clip length is approximately equal to aligned length
     * 
     * @param softClipLength total soft clipped bases
     * @param alignedLength total aligned bases
     * @return true if lengths are within tolerance
     */
    boolean lengthsMatch(int softClipLength, int alignedLength) {
        if (softClipLength == 0 || alignedLength == 0) {
            return false
        }
        
        double ratio = (double) softClipLength / (double) alignedLength
        return ratio >= (1.0 - lengthTolerance) && ratio <= (1.0 + lengthTolerance)
    }
    
    /**
     * Extract sequence from a read for the specified region
     * 
     * @param record SAMRecord
     * @param start start position (0-based, inclusive)
     * @param end end position (0-based, exclusive)
     * @return sequence as String
     */
    String extractSequence(SAMRecord record, int start, int end) {
        byte[] bases = record.getReadBases()
        if (start < 0 || end > bases.length || start >= end) {
            return ""
        }
        
        byte[] subBases = new byte[end - start]
        System.arraycopy(bases, start, subBases, 0, end - start)
        return new String(subBases)
    }
    
    /**
     * Align two sequences and return the alignment score as fraction of matching bases
     * 
     * @param seq1 first sequence
     * @param seq2 second sequence (will be reverse complemented)
     * @return fraction of matching bases (0.0 to 1.0)
     */
    double alignSequences(String seq1, String seq2) {

        return Align.global(gapOpenPenalty: 1, seq1, seq2).score / seq1.length()
    }
    
    /**
     * Detect if a read is a duplex read
     * 
     * @param record SAMRecord to check
     * @return true if the read appears to be a duplex read
     */
    boolean isDuplexRead(SAMRecord record) {
        // Skip unmapped reads
        if (record.getReadUnmappedFlag()) {
            return false
        }
        
        // Check if CIGAR is available
        if (record.getCigar() == null) {
            return false
        }
        
        // Extract soft clip information
        SoftClipInfo info = extractSoftClipInfo(record)
        
        // Must have soft clips
        if (!info.hasSoftClips()) {
            return false
        }
        
        // Check if either the leading or trailing soft clip individually matches the aligned length
        boolean leadingMatches = lengthsMatch(info.leadingSoftClip, info.alignedBases)
        boolean trailingMatches = lengthsMatch(info.trailingSoftClip, info.alignedBases)
        
        if (!leadingMatches && !trailingMatches) {
            return false
        }
        
        // Determine which end has the duplex soft clip (prefer the one that matches)
        boolean useTrailing
        if (trailingMatches && leadingMatches) {
            // Both match - prefer the larger one
            useTrailing = info.trailingSoftClip >= info.leadingSoftClip
        } else {
            useTrailing = trailingMatches
        }
        int softClipLength = useTrailing ? info.trailingSoftClip : info.leadingSoftClip
        
        // Extract sequences for comparison
        byte[] bases = record.getReadBases()
        if (bases == null || bases.length == 0) {
            return false
        }
        int readLength = bases.length
        
        String alignedSeq
        String softClipSeq
        
        if (useTrailing) {
            // Trailing soft clip: compare last N bp of aligned vs first N bp of soft clip
            int alignedEnd = readLength - softClipLength
            int alignedStart = Math.max(0, alignedEnd - alignmentWindowSize)
            alignedSeq = extractSequence(record, alignedStart, alignedEnd)
            
            int softClipStart = alignedEnd
            int softClipEnd = Math.min(readLength, softClipStart + alignmentWindowSize)
            softClipSeq = extractSequence(record, softClipStart, softClipEnd)
        } else {
            // Leading soft clip: compare last N bp of soft clip vs first N bp of aligned
            int softClipEnd = softClipLength
            int softClipStart = Math.max(0, softClipEnd - alignmentWindowSize)
            softClipSeq = extractSequence(record, softClipStart, softClipEnd)
            
            int alignedStart = softClipEnd
            int alignedEnd = Math.min(readLength, alignedStart + alignmentWindowSize)
            alignedSeq = extractSequence(record, alignedStart, alignedEnd)
        }
        
        // Align the sequences
        double alignmentScore = alignSequences(alignedSeq, softClipSeq)
        
        // Check if alignment exceeds threshold
        return alignmentScore >= alignmentThreshold
    }
}
