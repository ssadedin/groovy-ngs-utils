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
import htsjdk.samtools.BAMRecord
import htsjdk.samtools.SAMRecord

/**
 * A convenience class to streamline working with SAMRecord objects
 * 
 * @author simon.sadedin
 */
class SAMRecordCategory {
    @CompileStatic
    static List<XA> getAlternateAlignments(BAMRecord record) {
        String xas = record.getAttribute("XA")
        if(!xas)
            return null
        return xas.split(';').collect { String xa ->
            String [] parts = xa.split(',')
            new XA( chr: parts[0], pos: parts[1].substring(1) as Integer, cigar:parts[2],  nm : parts[-1] as Integer)
        }
    }
    
    @CompileStatic
    static Object asType(SAMRecord r, Class clazz) {
        new Region(r.referenceName, Math.min(r.alignmentStart, r.alignmentEnd)..Math.max(r.alignmentStart, r.alignmentEnd))
    }
    
    @CompileStatic
    static Object asType(BAMRecord r, Class clazz) {
        new Region(r.referenceName, Math.min(r.alignmentStart, r.alignmentEnd)..Math.max(r.alignmentStart, r.alignmentEnd))
    } 
    
    @CompileStatic
    static Object toRegion(BAMRecord r) {
        new Region(r.referenceName, Math.min(r.alignmentStart, r.alignmentEnd)..Math.max(r.alignmentStart, r.alignmentEnd))
    }  
}