package gngs;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.cram.encoding.readfeatures.BaseQualityScore;
import htsjdk.samtools.util.SequenceUtil;

public class CompactReadPair {
    
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
   
   public boolean isUnmapped() {
       return (r1ReferenceName == "*") || (r2ReferenceName == "*");
   }
   
   public void appendTo(StringBuilder b, SAMRecord r2, boolean addPosition) {
       
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
        b.append('@');
        b.append(r1Name);
        if(nameSuffix != null) { 
            b.append((CharSequence)nameSuffix);
        }        
        b.append(" 1:N:0:1\n");
        b.append(new String(basesAndQuals, 0, readLength));
        b.append("\n+\n");
        
        final int twoRL = readLength*2;
        for(int i=readLength; i<twoRL; ++i) {
            b.append(SAMUtils.phredToFastq(basesAndQuals[i]));
        } 
            
        // R2
        b.append("\n@");
        b.append(r1Name);
        if(nameSuffix != null) { 
            b.append((CharSequence)nameSuffix);
        }        
        
        b.append(" 2:N:0:1\n");
        b.append(SequenceUtil.reverseComplement(r2.getReadString()));
        b.append("\n+\n");
        
        byte [] bq = r2.getBaseQualities();
        for(int i=bq.length-1; i>=0; --i) {
            b.append(SAMUtils.phredToFastq(bq[i]));
        }
        b.append('\n');
   }
}
