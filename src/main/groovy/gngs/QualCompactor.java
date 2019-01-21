package gngs;

/**
 * A utility to compact quality scores by run length encoding them
 * <p>
 * Some specific optimisations are made to make things faster and more compact
 * for the case of base quality scores:
 * <li>We encode the length of the read as the first byte. To allow for values from
 *     0 - 256 we subtract 127 from the value since Java only allows signed values.
 *     Encoding the length of the read saves 
 * <li>The length of runs is similarly encoded as <code>value - 127</code>, which allows
 *     any given run of bases to be up to 256 long.
 * 
 * @author Simon Sadedin
 */
public class QualCompactor {
    
    /**
     * @param bases bases to compress
     * @return  byte array of bases, compacted to fit two bases per byte
     */
    final static byte [] compact(final byte [] bases) {
        final byte [] buffer = new byte[bases.length*2]; // if every base has different qual, we'll take twice the space
        int prev = bases[0];
        int count = -127;
        int bufferPos = 0;
        for(int i=0; i<bases.length; ++i) {
            final int b = bases[i];
            if(prev == b) {
                ++count;
            }
            else {
                buffer[bufferPos++] = (byte) count;
                buffer[bufferPos++] = (byte) prev;
                count = -126;
            }
            prev = b;
        }
        
        // We will be left with the last base unencoded
        if(count>-127) {
            buffer[bufferPos++] = (byte) count;
            buffer[bufferPos++] = (byte) prev;       
        }
                
        byte [] result = new byte[bufferPos+1];
        
        // Store the original length in the first 2 bytes
        result[0] = (byte)(bases.length - 127);
        
        System.arraycopy(buffer, 0, result, 1, bufferPos);
        return result;
    }
    
    final static byte [] expand(final byte [] compacted) {
        final int resultLength = compacted[0] + 127;
        byte [] result = new byte[resultLength];
        final int compactedLength = compacted.length;
        int resultPos = 0;
        for(int i=1; i<compactedLength; ++i) {
            int runLength = compacted[i++] + 127;
            byte value = compacted[i];
            for(int j=0; j<runLength;++j) {
                result[resultPos++] = value;
            }
        }
        return result;
    } 
}
