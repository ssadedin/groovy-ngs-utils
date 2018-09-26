/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2018 Simon Sadedin, ssadedin<at>gmail.com
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
package gngs.pair

import gngs.CompactReadPair
import gngs.ReadPair
import gngs.SAMRecordPair
import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import groovyx.gpars.actor.DefaultActor
import groovyx.gpars.util.DefaultMessageQueue
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SAMUtils
import htsjdk.samtools.util.SequenceUtil

/**
 * Converts reads to FASTQ format, and sends them in batches to a writer.
 * <p>
 * A {@link #buffer} is accumulated by appending each read requested to it in
 * FASTQ format, until the buffer exceeds a specified size. In general, this
 * buffer size should be optimised to ensure that the downstream I/O operations
 * are efficient, which generally means writing at least 10KB chunks at least,
 * but up to 1MB may be appropriate depending on the file system and OS.
 * 
 * @author Simon Sadedin
 */
@CompileStatic
@Log
class PairFormatter extends DefaultActor {
    
    StringBuilder buffer
    
    StringBuilder buffer2
    
    int maxBufferSize
    
    PairWriter writer
    
    PairWriter writer2
    
    int formatted = 0
    
    int buffered = 0
    
    String debugRead = null // "SN7001291:342:HFMC7BCXX:2:1113:1980:38527"
    
    /**
     * If true, the original position of the read will be appended 
     * to the read name. This allows backtracking of reads to compare between
     * alignments, etc.
     */
    boolean addPosition = false
    
    PairFormatter(int bufferSize, PairWriter writer, PairWriter writer2=null) {
        buffer = new StringBuilder(bufferSize+2000)
        this.maxBufferSize = bufferSize
        assert writer != null
        this.writer = writer
        this.writer2 = writer2
        if(writer2 != null)
            this.buffer2 = new StringBuilder(bufferSize+2000)
    }

    @CompileDynamic
    void act() {
        
        loop {
            react { Object msg ->
                assert msg != null
                if(msg == "stop") {
                    flushBuffer()
                    terminate()
                }
                else {
                    List msgList = (List)msg
                    process((ReadPair)msgList[0], (SAMRecord)msgList[1])
                }
                
            }
        }
    }
    
//        if(r2.getReadNegativeStrandFlag() && r2.getMateNegativeStrandFlag()) // read is aligned complemented already!
//        if(r2.getReadNegativeStrandFlag() && r2.getMateNegativeStrandFlag())  { // read is aligned reversed already! 
//            final int bqLength = bq.length;
//            for(int i=0; i<bqLength; ++i) {
//                b2.append(SAMUtils.phredToFastq(bq[i]));
//            }
//        }    
    
    void process(ReadPair pair, SAMRecord r2) {
        
        if(debugRead != null && (r2.readName == debugRead)) {
            log.info "Format: $debugRead"
        }
        
        if(pair instanceof CompactReadPair) {
            pair.appendTo(buffer, (StringBuilder) (buffer2 != null ? buffer2 : buffer), r2, this.addPosition)
        }
        else {
            ((SAMRecordPair)pair).appendTo(buffer, (buffer2 != null ? buffer2 : buffer), this.addPosition)
        }
        
        assert writer != null
        assert buffer != null
        
        if(buffer.size() > maxBufferSize) {
            flushBuffer()
        }
        else {
            buffered += 2
        }
    }
    
    void flushBuffer() {
        flushBufferAndWriter(buffer, writer)
        if(writer2 != null)
            flushBufferAndWriter(buffer2, writer2)
        formatted += buffered
        buffered = 0
    }
    
    void flushBufferAndWriter(StringBuilder buffer, PairWriter writer) {
        writer << [ content: buffer.toString(), reads: buffered ]
        writer.pending.addAndGet(buffered)
        buffer.setLength(0)
    }
}
