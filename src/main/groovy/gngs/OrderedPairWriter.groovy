/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2016 Simon Sadedin, ssadedin<at>gmail.com
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

package gngs

import groovy.transform.CompileStatic;
import htsjdk.samtools.SAMFileWriter

/**
 * When processing reads as pairs it can be tricky to write them out again while maintaining the
 * order expected for a "sorted" BAM file. This class buffers the second read until the next 
 * R1 is output that would violate the second read's sort order. 
 * 
 * <em>Note:</em> Currently this class filters out chimeric reads since ordering them is
 *                particularly tricky.
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class OrderedPairWriter implements Closeable {
    
    @Delegate
    SAMFileWriter samWriter 
    
    TreeSet<SAMRecordPair> buffer = new TreeSet({ SAMRecordPair p1, SAMRecordPair p2 ->
        p1.r2.alignmentStart.compareTo(p2.r2.alignmentStart)
    })
    
    int written = 0
    
    int currentReferenceIndex = -1
    
    ProgressCounter progress 
    
    OrderedPairWriter(SAMFileWriter w) {
        this.samWriter = w
        this.progress = new ProgressCounter(withTime:true, withRate:true, timeInterval:10000, extra: {
            "Written ${written}, Buffered ${buffer.size()} RefIndex=${currentReferenceIndex}"
        })
    }
    
    boolean verbose = false
    
    @CompileStatic
    void addAlignmentPair(SAMRecordPair pair) {
        
       if(pair.isChimeric())
           return
       
       progress.count()
       
       int r1ReferenceIndex = pair.r1.referenceIndex
       if(currentReferenceIndex < 0)     
           currentReferenceIndex = r1ReferenceIndex
       else
       if(r1ReferenceIndex != currentReferenceIndex) { // Switch between chromosomes, run down buffer
           flushBuffer()
       }
       
       currentReferenceIndex = r1ReferenceIndex
        
       SAMRecordPair bufferedPair = buffer.isEmpty() ? null : buffer.first() 
       while(bufferedPair != null && bufferedPair.r2.alignmentStart < pair.r1.alignmentStart) { 
           bufferedPair = buffer.pollFirst()
           if(verbose)
               println "WW: write R2 $bufferedPair.r1.readName ($bufferedPair.r1.referenceName:$bufferedPair.r1.alignmentStart,$bufferedPair.r2.referenceName:$bufferedPair.r2.alignmentStart)"
           samWriter.addAlignment(bufferedPair.r2) 
           ++written
           bufferedPair = buffer.isEmpty() ? null : buffer.first() 
       }
       
       if(verbose)
            println "WW: write R1 $pair.r1.readName ($pair.r1.referenceName:$pair.r1.alignmentStart,$pair.r2.alignmentStart)"
        
       samWriter.addAlignment(pair.r1)
       buffer << pair
    }
    
    @CompileStatic
    void flushBuffer() {
       for(SAMRecordPair bufferedPair in buffer) {
           samWriter.addAlignment(bufferedPair.r2)
           ++written
       }
       buffer.clear()
    }
    
    void close() {
        flushBuffer()
        this.samWriter.close()
        if(progress)
            progress.end()
    }
}
