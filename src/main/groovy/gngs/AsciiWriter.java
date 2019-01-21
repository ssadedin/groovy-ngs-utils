package gngs;

import java.io.IOException;
import java.io.OutputStream;
import java.io.Writer;

import htsjdk.samtools.util.StringUtil;

/**
 * 
 * Copied from HTSJDK but with constructor modified so that the buffer size can be customized
 * 
 * Fast (I hope) buffered Writer that converts char to byte merely by casting, rather than charset conversion.
 */
public final class AsciiWriter extends Writer {

    private final OutputStream os;
    // Buffer size has not been tuned.
    private final byte[] buffer;
    private int numBytes;

    /**
     * @param os need not be buffered as this class buffers
     */
    public AsciiWriter(final OutputStream os, int bufferSize) {
        this.os = os;
        numBytes = 0;
        this.buffer = new byte[bufferSize];
    }

    /**
     * flushes and closes underlying OutputStream.
     */
    @Override
    public void close() throws IOException {
        flush();
        os.close();
    }

    /**
     * flushes underlying OutputStream
     */
    @Override
    public void flush() throws IOException {
        os.write(buffer, 0, numBytes);
        numBytes = 0;
        os.flush();
    }

    /**
     * All other Writer methods vector through this, so this is the only one that must be overridden.
     */
    @Override
    public void write(final char[] chars, int offset, int length) throws IOException {
        while (length > 0) {
            final int charsToConvert = Math.min(length, buffer.length - numBytes);
            StringUtil.charsToBytes(chars, offset, charsToConvert, buffer, numBytes);
            numBytes += charsToConvert;
            offset += charsToConvert;
            length -= charsToConvert;
            if (numBytes == buffer.length) {
                os.write(buffer, 0, numBytes);
                numBytes = 0;
            }
        }
    }
}

