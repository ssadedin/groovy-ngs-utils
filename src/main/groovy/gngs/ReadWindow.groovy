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
package gngs

import groovy.transform.CompileStatic
import htsjdk.samtools.SAMRecord

/**
 * Models a window positioned centrally around a location in the genome containing overlapping reads
 * 
 * @author simon.sadedin
 */
@CompileStatic
class ReadWindow {
    
    /**
     * The of the window
     */
    int pos
    
    TreeMap<Integer, List<SAMRecord>> window = new TreeMap()
    
    String toString() {
        "Reads at ${window.count { it.value.size() > 1 }} positions: " + 
            window.grep { Map.Entry e -> ((List)e.value).size() > 1 }
    }
}
