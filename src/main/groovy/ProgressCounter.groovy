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

import groovy.transform.CompileStatic;

/**
 * Simple utility for displaying progress
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
                        System.err.println(new Date().toString() + ":\t Processed $count")
                    else
                        System.err.println "Processed $count"
                }
                lastPrintTimeMs = System.currentTimeMillis()
            }
        }
        ++count
    }
    
    void end() {
        System.err.println "Processed $count in ${System.currentTimeMillis() - startTimeMs} ms"
    }
}
