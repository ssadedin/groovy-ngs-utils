import java.util.Iterator;

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;

/**
 * Additional state about a read that is part of a pileup
 * 
 * @author simon.sadedin@mcri.edu.au
 */
public class PileupState {
    
    public SAMRecord read = null;
    
    public CigarElement cigarElement = null;
    
    public CigarOperator cigar = null;
    
    public Iterator<CigarElement> cigarIterator = null;
    
    /**
     * Offset within the current cigar element
     */
    public int cigarPos = -1;
    
    /**
     * Offset from start of read for the next 
     */
    public int nextReadPos = -1;
    
    public byte base;
    
    public byte qual;

    private int cigarLength;
    
    byte [] readBases = null;
    
    byte [] qualities = null;

    PileupState(SAMRecord r) {
        this.read = r;
        this.cigarIterator = r.getCigar().getCigarElements().listIterator();
        if(!cigarIterator.hasNext())
            throw new IllegalStateException("Bad CIGAR for read " + r.getReadName() + " has no elements");
        this.cigarElement = cigarIterator.next();
        this.cigarLength = cigarElement.getLength();
        this.cigar = cigarElement.getOperator();
        this.readBases = r.getReadBases();
        this.qualities = r.getBaseQualities();
        
        // Since soft clipped reads are before the alignment
        // start position, we have to skip over them
        // to be consistent with the overall pileup
        if(cigar == CigarOperator.S) {
            this.nextReadPos = cigarLength;
            this.cigarElement = cigarIterator.next();
            this.cigarLength = cigarElement.getLength();
            this.cigar = cigarElement.getOperator();
        }
        else {
          this.nextReadPos = 0;
        }
        this.cigarPos = -1;
    }
    
    /**
     * Move forward one reference position.
     * If the position covers a deletion, this does not alter the
     * read position. If the position covers an insertion, it will
     * jump forward over the insertion in the read.
     */
    public boolean next() {
        
        ++cigarPos;
        if(cigarPos >= cigarLength) {
            if(!cigarIterator.hasNext())
                return false;
            cigarElement = cigarIterator.next();
            cigarPos = 0;
            cigarLength = cigarElement.getLength();
            cigar = cigarElement.getOperator();
        }
        
        switch(cigar) {
            case D:
                base = 0;
                qual = 0;
                break;
                
            case I:
                
                // Skip the insertion
                nextReadPos += cigarElement.getLength();
                
                if(!cigarIterator.hasNext())
                    return false;
                
                base = readBases[nextReadPos];
                qual = qualities[nextReadPos];
                cigarElement = cigarIterator.next();
                cigar = cigarElement.getOperator();
                cigarLength = cigarElement.getLength();
                
                // This is set to zero because we are returning the
                // first element of the next cigar block
                cigarPos = 0;
                ++nextReadPos;
                break;
                
            default:
                if(nextReadPos >= readBases.length) {
                    System.err.println("Cigar exceeded length of read " + this.read.getReadName() + " in cigar op " + cigar + " at cigar position " + cigarPos);
                    base = 0;
                }
                else {
                  base = readBases[nextReadPos];
                  qual = qualities[nextReadPos];
                }
                ++nextReadPos;
        }
        
        return true;
    }
    
    /*
    public byte getBase() {
        switch(cigar) {
            case D:
                return -1;
                
            case I:
                throw new IllegalStateException("Pileup should never contain insertion");
                
            default:
                return readBases[readPos];
        }
    }
    */
}
