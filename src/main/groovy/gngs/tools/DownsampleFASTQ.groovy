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
package gngs.tools

import java.util.zip.GZIPOutputStream

import gngs.FASTQ
import gngs.ToolBase

/**
 * A simple tool for downsampling paired end FASTQ to achieve lower coverage.
 * 
 * @author Simon Sadedin
 */
class DownsampleFASTQ extends ToolBase {
    
    static void main(String[] args) {
        cli('DownsampleFASTQ -r <fraction of reads to preserve> -i1 <fastq r1> -i2 <fastq r2> -o1 <fastq output1> -o2 <fastq output2>', args) {
            i1 'FASTQ R1', args:1, required: true
            i2 'FASTQ R2', args:1, required: true
            o1 'Output FASTQ R2', args:1, required: true
            o2 'Output FASTQ R2', args:1, required: true
            r 'Fraction of reads to preserve', args:1, required: true
        } 
    }

    @Override
    public void run() {
        
        final double rate = opts.r.toDouble()
        Random random = new Random()
        OutputStream os1 = new GZIPOutputStream(new FileOutputStream(opts.o1))
        OutputStream os2 = new GZIPOutputStream(new FileOutputStream(opts.o2))
        os1.withWriter { w1 ->
            os2.withWriter { w2 ->
               FASTQ.filter(opts.i1, opts.i2, w1, w2) { r1, r2 ->
                    return (random.nextDouble()<rate)
                }            
            }
        }
     
    }

}
