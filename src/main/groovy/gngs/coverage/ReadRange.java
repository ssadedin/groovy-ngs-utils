package gngs.coverage;

import htsjdk.samtools.SAMRecord;

public class ReadRange {
    
    public ReadRange(SAMRecord r) {
        this.referenceIndex = r.getReferenceIndex();
        this.alignmentStart = r.getAlignmentStart();
        this.alignmentEnd = r.getAlignmentEnd();
        this.referenceName = r.getReferenceName();
    }
    
    final public int referenceIndex;
    
    final public int alignmentStart;
    
    final public int alignmentEnd;
    
    final public String referenceName;
}
