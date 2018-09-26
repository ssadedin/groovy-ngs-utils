package gngs;

import groovy.transform.CompileStatic;
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
   
   public void appendTo(StringBuilder r1Out, StringBuilder r2Out, SAMRecord r2, boolean addPosition) {
       
       // The logic here is complex
       // Read bases are all passed in relative to the forward strand of the reference
       // However we know that the two reads come from opposite strands, and hence
       // one read needs to be reverse complemented.
       //
       // We also know that the first read we encounter may be originally r1 or r2. 
       // So we also want to switch the order back so that R1 is actually written as 
       // the original R1, R2 as the original R2. This information comes from the 
       // firstOfPair / secondOfPair flags.
       // 
       // The canonical case is R1 -> <- R2
       // In that scenario, R1 is written as is, R2 is complemented
       // ie: (r1,r2')
       //
       // The next scenario is R2 -> <- R1
       // In this scenario we should switch the output buffer: since we encountered R2 first,
       // it will be stored as r1. However R1 will be have been complemented in alignment, 
       // as will have R2. So to get it right, we need to complement both. 
       // ie: (r2',r1')
       // 
       // Then there are more weird scenarios:
       // 
       // <-R1 R2->
       // This scenario is actually just the same as R1 -> <- R2, because both reads are 
       // actually still represented as bases wrt the forward strand. So this scenario just
       // comes back to:
       // (r1, r2')
       //
       // <-R2 <-R1
       
        final String r1Name = r2.getReadName();
        final int readLength = basesAndQuals.length >> 1;
        
//        StringBuilder b1 = r2.getFirstOfPairFlag() ? r1Out : r2Out;
//        StringBuilder b2 = r2.getFirstOfPairFlag() ? r2Out : r1Out;
        
        StringBuilder b1 = r1Out;
        StringBuilder b2 = r2Out;
                
        StringBuilder nameSuffix = null;
        
        if(addPosition) { 
            nameSuffix = new StringBuilder();
            nameSuffix.append(':');
            nameSuffix.append(r1ReferenceName);
            nameSuffix.append(':');
            nameSuffix.append(r1AlignmentStart);
            nameSuffix.append(':');
            if(isChimeric()) {
                nameSuffix.append(r2ReferenceName);
                nameSuffix.append(':');
            }
            nameSuffix.append(r2AlignmentStart);
        }
           
        //////////////////////////////////////////////////////////////////////
        // R1
        //////////////////////////////////////////////////////////////////////
        
        // r1 could be R1 or R2
        // If it's R1 and aligned to neg strand, then we shoudl flip back
        // if it's R2 and aligned to pos strand, then we sh
        final boolean flipR1 = r2.getMateNegativeStrandFlag();
        
        final byte [] basesAndQuals = this.basesAndQuals;
        b1.append('@');
        b1.append(r1Name);
        if(nameSuffix != null) { 
            b1.append((CharSequence)nameSuffix);
        }        
        b1.append(" 1:N:0:1\n");
        if(flipR1)
            SequenceUtil.reverseComplement(basesAndQuals, 0, readLength);
        b1.append(new String(basesAndQuals, 0, readLength));
        b1.append("\n+\n");
        
        final int twoRL = readLength*2;
        
        if(flipR1) { // read is aligned complemented, we have to turn it back
            for(int i=twoRL-1; i>=readLength; --i) {
                b1.append(SAMUtils.phredToFastq(basesAndQuals[i]));
            }  
        }
        else
        for(int i=readLength; i<twoRL; ++i) {
            b1.append(SAMUtils.phredToFastq(basesAndQuals[i]));
        } 
        b1.append('\n');
            
        //////////////////////////////////////////////////////////////////////
        // R2
        //////////////////////////////////////////////////////////////////////
        b2.append("@");
        b2.append(r1Name);
        if(nameSuffix != null) { 
            b2.append((CharSequence)nameSuffix);
        }        
        
        b2.append(" 2:N:0:1\n");
        
        final boolean flipR2 = r2.getReadNegativeStrandFlag();
        if(flipR2)
            b2.append(SequenceUtil.reverseComplement(r2.getReadString()));
        else
            b2.append(r2.getReadString());
        
        b2.append("\n+\n");
        
        byte [] bq = r2.getBaseQualities();
        if(flipR2) {
            for(int i=bq.length-1; i>=0; --i) {
                b2.append(SAMUtils.phredToFastq(bq[i]));
            }
        }
        else {
            for(int i=0; i<bq.length; ++i) {
                b2.append(SAMUtils.phredToFastq(bq[i]));
            }            
        }
        b2.append('\n');
   }
}
