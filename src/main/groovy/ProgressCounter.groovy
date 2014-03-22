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
    
    int lineInterval = 500
    
    long timeInterval = 30000
    
    long lastPrintTimeMs = System.currentTimeMillis()
    
    long startTimeMs = -1L
    
    boolean withTime = false
    
    String prefix = null
    
    PrintStream out = System.err
    
    public ProgressCounter() {
    }
    
    @CompileStatic
    void count(Closure c = null) {
        if(startTimeMs < 0)
            startTimeMs = System.currentTimeMillis()
        if(count % lineInterval == 0) {
            if(System.currentTimeMillis() - lastPrintTimeMs > timeInterval) {
                if(c!=null)
                  c(count)
                else {
                    if(withTime)
                        out.println(new Date().toString() + (prefix?"\t$prefix":"") + " :\t Processed $count")
                    else
                        out.println((prefix?"\t$prefix :\t":"") + "Processed $count")
                }
                lastPrintTimeMs = System.currentTimeMillis()
            }
        }
        ++count
    }
    
    void end() {
        if(count == 0) {
            out.println "0 records processed"
            return
        }
            
        long endTime = System.currentTimeMillis()
        def timeDelta = TimeCategory.minus(new Date(endTime),new Date(startTimeMs))
        out.println "Processed $count in ${timeDelta} @ ${1000*count/((endTime-startTimeMs)+1)} per second"
    }
    
    void getAbort() { 
        throw new Abort() 
    }; 

    static withProgress(Closure c) { 
        ProgressCounter progress = new ProgressCounter(withTime:true)
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
