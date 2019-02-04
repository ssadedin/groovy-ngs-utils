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
import groovy.transform.Memoized
import groovy.transform.ToString
import htsjdk.samtools.Cigar
import htsjdk.samtools.TextCigarCodec

@CompileStatic
@ToString(includeNames=true)
class AltContig {
    
    AltContig(String name, String chr, int start, String cigarString) {
        this.name = name
        this.cigarString = cigarString
        
        
        Cigar cigar = TextCigarCodec.decode(cigarString)
        
        int cigarSize = cigar.getReferenceLength()
        this.primaryRegion = new Region(chr, start..(start + cigarSize), name: name, cigar: cigarString)
    }
    
    String name
    
    String cigarString
    
    Region primaryRegion
}

/**
 * Parses an alt-index in SAM format as provided with GRCh38 to define the alt contigs,
 * including extracting the regions in the primary assembly covered by each alt contig
 * <p>
 * Example Usage:
 * <pre>
 * println "chr17_KI270908v1_alt overlaps " + new AltContigIndex("Homo_sapiens_assembly38.fasta.64.alt").contigs.chr17_KI270908v1_alt.primaryRegion + " in the GRCh38 reference"
 * </pre>
 * 
 * @author Simon Sadedin
 */
class AltContigIndex {
    
    HashMap<String,AltContig> contigs = null
    
    AltContigIndex(fileLike) {
        this.contigs = Utils.reader(fileLike).withReader { r ->
            r.readLines()*.tokenize().grep {!it[0].startsWith('@SQ') &&  it[2] != '*' }collectEntries { List parts ->
                [parts[0], new AltContig(parts[0], parts[2], parts[3]?.toInteger(), parts[5])]
            }
        }
    }
    
    @Memoized
    Regions getRegions() {
        this.contigs*.value*.primaryRegion as Regions
    }
}
