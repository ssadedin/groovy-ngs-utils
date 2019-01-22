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

import java.io.Writer
import java.util.Map
import java.util.concurrent.atomic.AtomicInteger
import gngs.RegulatingActor
import gngs.SAMRecordPair
import groovy.transform.CompileStatic
import groovyx.gpars.actor.DefaultActor

/**
 * A simple class that whose purpose is purely to paralleise
 * writing output to a Writer. 
 * <p>
 * Each message contains two attributes:
 * <li>content - a String value representing the data to be written
 * <li>reads   - the count of reads contained in the data
 * Each message can contain many reads and therefore potentially represents a large 
 * amount of data. Therfore the queuing on this queue should be constrained to be quite
 * small. For example, if the buffer size is 1MB then 50 queued messages will represent
 * 50MB of memory.
 * <p>
 * A {@link #pending} count is decremented with each batch of reads written.
 * This count needs to be incremented by the caller for it to be accurate. If 
 * this is done, the pending count will reflect how many reads are queued, waiting 
 * to be written, and thus give an estimate of how much memory is currently 
 * being consumed by buffered reads.
 * 
 * @author Simon Sadedin
 */
@CompileStatic
class PairWriter extends RegulatingActor<Map<String,Object>> {
    
    Writer out
    
    int written = 0
    
    public AtomicInteger pending = new AtomicInteger()
    
    PairWriter(Writer writer, int queueSize) {
        super(queueSize,queueSize*2)
        this.out = writer
    }
    
    @Override
    public void process(Map<String, Object> msg) {
        String value = (String)msg['content']
        out.write(value, 0, value.size())
        int readCount = (int)msg['reads']
        written = written + readCount
        pending.addAndGet(-readCount)
    }
}


