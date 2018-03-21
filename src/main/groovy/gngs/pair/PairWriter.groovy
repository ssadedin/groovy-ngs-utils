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
import java.util.concurrent.atomic.AtomicInteger
import gngs.SAMRecordPair
import groovy.transform.CompileStatic
import groovyx.gpars.actor.DefaultActor

/**
 * A simple class that whose purpose is purely to paralleise
 * writing output to a Writer.
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
class PairWriter extends DefaultActor {
    
    Writer out
    
    int written = 0
    
    public AtomicInteger pending = new AtomicInteger()
    
    PairWriter(Writer writer) {
        this.out = writer
    }
    
    void act() {
        loop {
            react { msg ->
                if(msg == "stop")
                    terminate()
                else {
                    Map msgMap = (Map)msg
                    out.append((String)msg['content'])
                    int readCount = (int)msg['reads']
                    written = written + readCount
                    pending.addAndGet(-readCount)
                }
            }
        }
    }
}


