/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2013 Simon Sadedin, ssadedin<at>gmail.com
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

import groovy.time.TimeCategory;
import groovy.transform.CompileStatic;
import java.util.logging.Logger

/**
 * Simple utility for displaying progress 
 * <p>
 * The code doing the 'work' must call 'count' on the progress counter.
 * This will print out a progress report at spaced repetitions (the
 * default is every 500 counts or 30 seconds, whichever happens first).
 * <p>
 * Optionally, the 'work' code can call "end()" when it is finished,
 * which will then print out some brief summary statistics for the whole
 * operation.
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class ProgressCounter {
    
    int count = 0
    
    int lastPrintCount = 0
    
    int lineInterval = 500
    
    long timeInterval = 30000
    
    long lastPrintTimeMs = System.currentTimeMillis()
    
    long startTimeMs = -1L
    
    boolean withTime = false
    
    boolean withRate = false
    
    String prefix = null
    
    PrintStream out = System.err
    
    Logger log
    
    Closure extra = null
    
    @CompileStatic
    void count(Closure c = null) {
        if(startTimeMs < 0)
            startTimeMs = System.currentTimeMillis()
        if(count % lineInterval == 0) {
            long deltaMs = System.currentTimeMillis() - lastPrintTimeMs
            if(deltaMs > timeInterval) {
                if(c!=null)
                  c(count)
                else {
                    
                    String rateInfo = ""
                    if(withRate) {
                        float rate = 1000.0*((count-lastPrintCount)/((float)(deltaMs+1)))
                        rateInfo = String.format(" @ %.2f/s ", rate)
                    }
                        
                   String extraInfo
                   if(extra != null) {
                       extraInfo = extra()
                   }   
                    
                   String progressLine
                   if(withTime)
                       progressLine = new Date().toString() + (prefix?"\t$prefix":"") + " :\t Processed $count" + rateInfo + (extraInfo? " " + extraInfo : "")
                   else
                       progressLine = (prefix?"\t$prefix :\t":"") + "Processed ${String.valueOf(count).padRight(12)}" + rateInfo + (extraInfo? " " + extraInfo : "")
                       
                       
                   if(log != null)
                       log.info(progressLine)
                   else
                       out.println(progressLine) 
                }
                lastPrintTimeMs = System.currentTimeMillis()
                lastPrintCount = count
            }
        }
        ++count
    }
    
    void end() {
        if(count == 0) {
            if(log)
                log.info "0 records processed"
            else
                out.println "0 records processed"
            return
        }
        
        String extraInfo = ""
        if(extra != null) {
            extraInfo = extra()
        }
        
        long endTime = System.currentTimeMillis()
        long deltaMs = (endTime-startTimeMs)+1 
        def timeDelta = TimeCategory.minus(new Date(endTime),new Date(startTimeMs))
        float rate = 1000.0d*((count)/((float)(deltaMs+1)))
        String rateInfo = String.format(" @ %.2f/s ", rate)
        
        String progressLine = "Processed ${String.valueOf(count)} in ${timeDelta} ${rateInfo} " + (extra != null ? " ($extraInfo)" : "")
        if(log)
            log.info(progressLine)
        else
            out.println(progressLine)
    }
    
    void getAbort() { 
        throw new Abort() 
    }; 

    static withProgress(Map options=[:],Closure c) { 
        ProgressCounter progress = new ProgressCounter()
        
        for(option in options) {
            progress[option.key] = option.value
        }
        
        c.delegate = progress
        try {
            c(progress)
        }
        catch(Abort e) {
            // Ignore
        }
        finally {
          progress.end()  
        }
    }
}
