package gngs;

/**
 * A utility to compact bases to a two-per-byte representation
 * 
 * @author Simon Sadedin
 */
public class BaseCompactor {
    
    public static final byte T = (byte)'T';
    public static final byte A = (byte)'A';
    public static final byte C = (byte)'C';
    public static final byte G = (byte)'G';
    public static final byte N = (byte)'N';
       
            
    private final static byte [] COMPACT_UNMAPPING = new byte[] {
        0, // Filler: we reserve zero as the flag to detect unitialized / empty slots
        A,
        T,
        C,
        G,
        N
    };
    
    private final static byte [] COMPACT_MAPPING = new byte[85];
    
    static {
        COMPACT_MAPPING[A] = 1;
        COMPACT_MAPPING[T] = 2;
        COMPACT_MAPPING[C] = 3;
        COMPACT_MAPPING[G] = 4;
        COMPACT_MAPPING[N] = 5;
    }
    
    /**
     * Bases are packed two-per-byte resulting in effectively a 4-bit representation.
     * <p>
     * The last byte could have either 2 bases or 1 base in it. This is recognisable
     * because the zero value is not a valid base, so if the high 4 bits of the last byte
     * are zero, then there was an odd number of bases to store.
     * 
     * @param bases bases to compress
     * @return  byte array of bases, compacted to fit two bases per byte
     */
    final static byte [] compact(final byte [] bases) {
        
        final int blen = (bases.length>>1)<<1; // round down to even
        final byte [] result = new byte[(bases.length>>1)+(bases.length%2)];
        int resultIndex = 0;
        for(int i=0; i<blen; ++i) {
           result[resultIndex] = (byte)(COMPACT_MAPPING[bases[i]] | (COMPACT_MAPPING[bases[++i]]<<4));
           ++resultIndex;
        }
        if(bases.length>blen) {
            result[resultIndex] = COMPACT_MAPPING[bases[bases.length-1]]; // odd number of bases
        }
            
        return result;
    }
    
    final static byte [] expand(final byte [] compacted) {
        final int compactedLength = compacted.length;
        final int safeCompactedLength = compactedLength-1; // there is special logic for the last byte, so don't do that in the loop
        final byte lastByte = compacted[safeCompactedLength];
        final int resultLength = (lastByte & 0xF0)>0 ? (compactedLength * 2) : (safeCompactedLength * 2 + 1);
        final byte [] result = new byte [resultLength];
        
        int resultIndex = 0;
        for(int i=0; i<safeCompactedLength; ++i) {
            final byte encoded = compacted[i];
            result[resultIndex++] = COMPACT_UNMAPPING[encoded & 0x0F];
            result[resultIndex++] = COMPACT_UNMAPPING[encoded >> 4];
        }
        result[resultIndex++] = COMPACT_UNMAPPING[lastByte & 0x0F];
        if((lastByte & 0xF0) > 0) {
            result[resultIndex] = COMPACT_UNMAPPING[(lastByte & 0xF0) >> 4];
        }
        
        return result;
    } 
}
