package gngs;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.cram.encoding.readfeatures.BaseQualityScore;
import htsjdk.samtools.util.SequenceUtil;

public class CompactReadPair implements ReadPair {
    
    public String r1ReferenceName;
    
    public String r2ReferenceName;
    
    public int r1AlignmentStart;
    
    public int r2AlignmentStart;
    
    public byte [] basesAndQuals;
    
    public CompactReadPair(SAMRecord read) {
        // Note that picard either uses a common string from the sequence dictionary for
        // all reads OR interns the string, so there isn't any point trying that here
        r1ReferenceName = read.getReferenceName();
        r2ReferenceName = read.getMateReferenceName();
        if(r1ReferenceName.equals(r2ReferenceName))
            r2ReferenceName = null;
        r1AlignmentStart = read.getAlignmentStart();
        r2AlignmentStart = read.getMateAlignmentStart();
        
        int readLength = read.getReadLength();
        basesAndQuals = new byte[readLength*2];
        
        System.arraycopy(read.getReadBases(),0,basesAndQuals,0,readLength);
        System.arraycopy(read.getBaseQualities(),0,basesAndQuals,readLength,readLength);
   }
    
   public int getReadLength() {
       return basesAndQuals.length >> 1;
   }
    
   public boolean isChimeric() {
       return r2ReferenceName != null;
   }
   
   public boolean getUnmapped() {
       return (r1ReferenceName == "*") || (r2ReferenceName == "*");
   }
   
   public boolean notInRegions(Regions regions) {
       final String chr1 = this.r1ReferenceName;
       final int r1p = this.r1AlignmentStart;
       final int len = this.getReadLength();
       if(regions.overlaps(chr1, r1p , r1p+len))
           return false;
       final String chr2 = this.r2ReferenceName != null ? this.r2ReferenceName:chr1;
       final int r2p = this.r2AlignmentStart;
       if(regions.overlaps(chr2, r2p , r2p+len)) {
           return false;
       }
       return true;
   }
   
   public void appendTo(StringBuilder b1, StringBuilder b2, SAMRecord r2, boolean addPosition) {
       
        String r1Name = r2.getReadName();
        int readLength = basesAndQuals.length >> 1;
        
        StringBuilder nameSuffix = null;
        
        if(addPosition) { 
            nameSuffix = new StringBuilder();
            nameSuffix.append(':');
            nameSuffix.append(r1ReferenceName);
            nameSuffix.append(':');
            nameSuffix.append(r1AlignmentStart);
            nameSuffix.append(':');
            if(isChimeric())
                nameSuffix.append(r2ReferenceName);
            nameSuffix.append(r2AlignmentStart);
        }
           
        // R1
        final byte [] basesAndQuals = this.basesAndQuals;
        b1.append('@');
        b1.append(r1Name);
        if(nameSuffix != null) { 
            b1.append((CharSequence)nameSuffix);
        }        
        b1.append(" 1:N:0:1\n");
        b1.append(new String(basesAndQuals, 0, readLength));
        b1.append("\n+\n");
        
        final int twoRL = readLength*2;
        for(int i=readLength; i<twoRL; ++i) {
            b1.append(SAMUtils.phredToFastq(basesAndQuals[i]));
        } 
        b1.append('\n');
            
        // R2
        b2.append("@");
        b2.append(r1Name);
        if(nameSuffix != null) { 
            b2.append((CharSequence)nameSuffix);
        }        
        
        b2.append(" 2:N:0:1\n");
        b2.append(SequenceUtil.reverseComplement(r2.getReadString()));
        b2.append("\n+\n");
        
        byte [] bq = r2.getBaseQualities();
        for(int i=bq.length-1; i>=0; --i) {
            b2.append(SAMUtils.phredToFastq(bq[i]));
        }
        b2.append('\n');
   }
}
