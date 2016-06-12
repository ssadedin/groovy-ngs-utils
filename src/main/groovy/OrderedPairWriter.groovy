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
    
    OrderedPairWriter(SAMFileWriter w) {
        this.samWriter = w
    }
    
    boolean verbose = false
    
    @CompileStatic
    void addAlignmentPair(SAMRecordPair pair) {
        
       if(verbose)
            println "WW: write R1 ($pair.r1.referenceName:$pair.r1.alignmentStart,$pair.r2.alignmentStart)"
        
       if(pair.isChimeric())
           return
        
       SAMRecordPair bufferedPair = buffer.isEmpty() ? null : buffer.first() 
       while(bufferedPair != null && bufferedPair.r2.alignmentStart < pair.r1.alignmentStart) {
           bufferedPair = buffer.pollFirst()
           if(verbose)
               println "WW: write R2 ($pair.r1.referenceName:$pair.r1.alignmentStart,$pair.r2.referenceName:$pair.r2.alignmentStart)"
           samWriter.addAlignment(bufferedPair.r2) 
           bufferedPair = buffer.isEmpty() ? null : buffer.first() 
       }
       
       samWriter.addAlignment(pair.r1)
       buffer << pair
    }
    
    void close() {
        for(SAMRecordPair pair in buffer) {
            this.samWriter.addAlignment(pair.r2)
        }
        this.samWriter.close()
    }
}
