/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) Simon Sadedin, ssadedin<at>gmail.com
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
package gngs;

import java.io.IOException;
import java.util.concurrent.atomic.AtomicLong;
import java.util.logging.Logger;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.xerial.snappy.Snappy;

import gngs.pair.Paired;
import graxxia.IntegerStats;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.SequenceUtil;

/**
 * Stores one read of a pair and links to the second read in a way that offers
 * a useful amount of compression for applications that need to store many
 * reads in memory.
 * <p>
 * By default, reads are compressed using a combination of two compression methods:
 * <p>
 * <li>The snappy compression for base qualities
 * <li>Packing of bases to two-per-byte for bases
 * <p>
 * Alternative methods for compression are configurable, however for most use cases, these
 * methods empirically outperform others in my experience. In particular, snappy compression 
 * of raw bases does not outperform simply packing them two-per-byte.
 * <p>
 * Note: 2-bit (4-per-byte) representation is not possible because we have to encode N bases.
 * Theoretically, if alignment to byte boundaries is sacrificed then 3-bit representation could
 * be achieved, but this has not been tested.
 * 
 * @author Simon Sadedin
 */
public class CompactReadPair implements ReadPair {
    
    private static Logger log = Logger.getLogger("CompactReadPair");
    
    public String r1ReferenceName;
    
    public String r2ReferenceName;
    
    public int r1AlignmentStart;
    
    public int r2AlignmentStart;
    
    public byte [][] compressedBases;
    
    private final short readLength;
    
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
    
    public static IntegerStats baseCompressionStats = new IntegerStats(200);
    public static IntegerStats qualCompressionStats = new IntegerStats(200);
    
    static ReadCompressor initCompressor() {
        String compressorType = System.getProperty("bazam.compressor", "snappyhybrid");
        log.info("Reads will be compressed using compression method: " + compressorType);
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
    
    /**
     * Create a compact read pair that encodes the default base qualities
     * <p>
     * Note: the first read is initialised by this constructor. The second read is not stored itself, rather details
     * are extracted from the metadata of the first read about its position. Usage of this class as a "pair" requires
     * provision of the second read explicitly.
     * 
     * @param read  the first read to encode
     * @throws IOException
     */
    public CompactReadPair(final SAMRecord read) throws IOException {
        this(read, null);
    }
    
    /**
     * Create a compact read pair that encodes base qualities extracted from the given tag.
     * 
     * @param read
     * @param baseQualityTag
     * @throws IOException
     */
    public CompactReadPair(final SAMRecord read, final String baseQualityTag) throws IOException {
        // Note that picard either uses a common string from the sequence dictionary for
        // all reads OR interns the string, so there isn't any point trying that here
        r1ReferenceName = read.getReferenceName();
        r2ReferenceName = read.getMateReferenceName();
        if(r1ReferenceName.equals(r2ReferenceName))
            r2ReferenceName = null;
        
        r1AlignmentStart = read.getAlignmentStart();
        r2AlignmentStart = read.getMateAlignmentStart();
        
        readLength = (short)read.getReadLength();
       
        final byte[] quals = (baseQualityTag != null) ? SAMUtils.fastqToPhred(read.getStringAttribute(baseQualityTag)) : read.getBaseQualities();
        
        final byte[] bases = read.getReadBases();
        compressedBases = compressor.compress(bases, quals);
        
        long mem = compressedBases[0].length;
        if(compressedBases.length>1)
            mem += compressedBases[1].length;
        
        
        memoryStats.addValue(memoryUsage.addAndGet(mem));
        
        baseCompressionStats.addIntValue((int)(100d * compressedBases[0].length / (double)readLength));
        if(compressedBases.length>1)
            qualCompressionStats.addIntValue((int)(100d * compressedBases[1].length / (double)readLength));
        
        currentCount.incrementAndGet();
   }
    
   public int getReadLength() {
       return (int)readLength;
   }
    
   public boolean isChimeric() {
       return r2ReferenceName != null;
   }
   
   public boolean getUnmapped() {
       return (r1ReferenceName.equals("*") || r2ReferenceName.equals("*"));
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
   
   public void appendTo(final String r1Name, final StringBuilder r1Out, final StringBuilder r2Out, Paired pairInfo, boolean addPosition) throws IOException {
       
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
        
        CompactReadPair r2 = (CompactReadPair)pairInfo.r2;
        
        // r1 could be R1 or R2
        // If it's R1 and aligned to neg strand, then we should flip back
        // if it's R2 and aligned to pos strand, then we sh
        final boolean flipR1 = pairInfo.getR1NegativeStrandFlag(); // TODO: FIXME r2.getMateNegativeStrandFlag();
       
        byte [][] expanded = compressor.expand(this.compressedBases);
        byte [] decompressedBases = expanded[0];
        byte [] decompressedQuals = expanded[1];
        if(flipR1)
            SequenceUtil.reverseComplement(decompressedBases);
         
        boolean reordered = pairInfo.getReorderRequired();
        if(reordered) {
            StringBuilder tmp = b1;
            b1 = b2;
            b2 = tmp;
            
            b2.append("@");
            b2.append(r1Name);
            if(nameSuffix != null) { 
                b2.append((CharSequence)nameSuffix);
            }        
            appendR2(pairInfo, b2, r2);
        }

        b1.append('@');
        b1.append(r1Name);
        if(nameSuffix != null) { 
            b1.append((CharSequence)nameSuffix);
        }        
        b1.append(" 1:N:0:1\n");
       
        b1.append(new String(decompressedBases, 0, readLength));
        b1.append("\n+\n");
        
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
       
        if(!reordered) {
            b2.append("@");
            b2.append(r1Name);
            if(nameSuffix != null) { 
                b2.append((CharSequence)nameSuffix);
            }        
            appendR2(pairInfo, b2, r2);
        }
        
        long mem = -compressedBases[0].length;
        mem -= r2.compressedBases[0].length;

        if(compressedBases.length>1) {
            mem -= compressedBases[1].length;
            mem -= r2.compressedBases[1].length;
        }

        memoryUsage.addAndGet(mem);
        currentCount.decrementAndGet();
   }

	private void appendR2(final Paired pairInfo, final StringBuilder b2, final CompactReadPair r2) throws IOException {
		b2.append(" 2:N:0:1\n");
		final byte [][] expandedR2 = compressor.expand(r2.compressedBases);
		final byte [] decompressedR2Bases = expandedR2[0];
		final byte [] decompressedR2Quals = expandedR2[1];
		final boolean flipR2 = pairInfo.getR2NegativeStrandFlag(); 
		if(flipR2)
			SequenceUtil.reverseComplement(decompressedR2Bases);
		b2.append(new String(decompressedR2Bases,0,decompressedR2Bases.length));
		
		b2.append("\n+\n");
		
		if(flipR2) {
			for(int i=decompressedR2Quals.length-1; i>=0; --i) {
				b2.append(SAMUtils.phredToFastq(decompressedR2Quals[i]));
			}
		}
		else {
			for(int i=0; i<decompressedR2Quals.length; ++i) {
				b2.append(SAMUtils.phredToFastq(decompressedR2Quals[i]));
			}            
		}
		b2.append('\n');
	}
}
