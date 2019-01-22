package gngs;

import java.io.IOException;
import java.util.concurrent.atomic.AtomicLong;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.xerial.snappy.Snappy;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.SequenceUtil;

public class CompactReadPair implements ReadPair {
    
    public String r1ReferenceName;
    
    public String r2ReferenceName;
    
    public int r1AlignmentStart;
    
    public int r2AlignmentStart;
    
    public byte [][] compressedBases;
    
    short readLength;
    
    public static interface ReadCompressor {
        byte [][] compress(byte[] bases, byte[] quals) throws IOException;
        byte [][] expand(byte[][] basesAndQuals) throws IOException;
    }
    
    
    static class SnappyReadCompressor implements ReadCompressor {

        @Override
        public byte[][] compress(final byte[] bases, final byte [] quals) throws IOException {
            return new byte[][] { Snappy.compress(bases), Snappy.compress(quals) };
        }

        @Override
        public byte[][] expand(byte[][] basesAndQuals) throws IOException {
            return new byte[][]{Snappy.uncompress(basesAndQuals[0]), Snappy.uncompress(basesAndQuals[1])};
        }
    }
    
    final static class SnappyCompactorHybridCompressor implements ReadCompressor {

        @Override
        public byte[][] compress(byte[] bases, byte[] quals) throws IOException {
            return new byte[][] { BaseCompactor.compact(bases), Snappy.compress(quals) };
        }

        @Override
        public byte[][] expand(byte[][] basesAndQuals) throws IOException {
            return new byte[][]{BaseCompactor.expand(basesAndQuals[0]), Snappy.uncompress(basesAndQuals[1])};
        }
    } 
    
    
    final static class CompactorCompressor implements ReadCompressor {

        @Override
        public byte[][] compress(final byte[] bases, final byte[] quals) throws IOException {
            return new byte[][] { BaseCompactor.compact(bases), QualCompactor.compact(quals) };
        }

        @Override
        public byte[][] expand(final byte[][] basesAndQuals) throws IOException {
            return new byte[][]{BaseCompactor.expand(basesAndQuals[0]), QualCompactor.expand(basesAndQuals[1])};
        }
    } 
    
    final static class SharedArrayCompressor implements ReadCompressor {

        @Override
        public byte[][] compress(final byte[] bases, final byte[] quals) throws IOException {
            final int readLength = (short)bases.length;
            byte [] basesAndQuals = new byte[readLength*2];
            System.arraycopy(bases,0,basesAndQuals,0,readLength);
            System.arraycopy(quals,0,basesAndQuals,readLength,readLength);
            return new byte[][] { basesAndQuals };
        }

        @Override
        public byte[][] expand(final byte[][] basesAndQuals) throws IOException {
            final byte [] data = basesAndQuals[0];
            final int readLength = data.length>>1;
            final byte[][] results = new byte[][] {
                new byte[readLength],
                new byte[readLength]
            };
            
            System.arraycopy(data, 0, results[0], 0, readLength);
            System.arraycopy(data, readLength, results[1], 0, readLength);            
            
            return results;
        }
    }    
    
    static AtomicLong memoryUsage = new AtomicLong(0);
    
    static AtomicLong currentCount = new AtomicLong(0);
    
    public static SummaryStatistics memoryStats  = new SummaryStatistics();
    
    static ReadCompressor initCompressor() {
        String compressorType = System.getProperty("bazam.compressor", "snappyhybrid");
        if(compressorType.equals("snappyhybrid")) {
            return new SnappyCompactorHybridCompressor();
        }
        else
        if(compressorType.equals("gngscompactor")) {
            return new CompactorCompressor();
        }
        else
        if(compressorType.equals("snappy")) {
            return new SnappyReadCompressor();
        }
        else
        if(compressorType.equals("none")) {        
            return new SharedArrayCompressor();
        }
        else {
            throw new IllegalStateException("bazam.compressor property was set to invalid value " + compressorType + " please choose a valid type, eg: gngscompactor");
        }
    }
    
    final static ReadCompressor compressor = initCompressor();
    
    public CompactReadPair(SAMRecord read) throws IOException {
        // Note that picard either uses a common string from the sequence dictionary for
        // all reads OR interns the string, so there isn't any point trying that here
        r1ReferenceName = read.getReferenceName();
        r2ReferenceName = read.getMateReferenceName();
        if(r1ReferenceName.equals(r2ReferenceName))
            r2ReferenceName = null;
        r1AlignmentStart = read.getAlignmentStart();
        r2AlignmentStart = read.getMateAlignmentStart();
        
        readLength = (short)read.getReadLength();
       
        compressedBases = compressor.compress(read.getReadBases(), read.getBaseQualities());
        
        long mem = compressedBases[0].length;
        if(compressedBases.length>1)
            mem += compressedBases[1].length;
        
        memoryStats.addValue(memoryUsage.addAndGet(mem));
        
        currentCount.incrementAndGet();
   }
    
   public int getReadLength() {
       return (int)readLength;
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
   
   public void appendTo(StringBuilder r1Out, StringBuilder r2Out, SAMRecord r2, boolean addPosition) throws IOException {
       
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
        
        byte [][] expanded = compressor.expand(this.compressedBases);
        byte [] decompressedBases = expanded[0];
        byte [] decompressedQuals = expanded[1];
        
        b1.append('@');
        b1.append(r1Name);
        if(nameSuffix != null) { 
            b1.append((CharSequence)nameSuffix);
        }        
        b1.append(" 1:N:0:1\n");
        if(flipR1)
            SequenceUtil.reverseComplement(decompressedBases);
        
        b1.append(new String(decompressedBases, 0, readLength));
        b1.append("\n+\n");
        
        final int twoRL = readLength*2;
        
        if(flipR1) { // read is aligned complemented, we have to turn it back
            for(int i=readLength-1; i>=0; --i) {
                b1.append(SAMUtils.phredToFastq(decompressedQuals[i]));
            }  
        }
        else
        for(int i=0; i<readLength; ++i) {
            b1.append(SAMUtils.phredToFastq(decompressedQuals[i]));
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
        
        long mem = -compressedBases[0].length;
        if(compressedBases.length>1)
            mem -= compressedBases[1].length;
        
        memoryUsage.addAndGet(mem);
        currentCount.decrementAndGet();
   }
}
