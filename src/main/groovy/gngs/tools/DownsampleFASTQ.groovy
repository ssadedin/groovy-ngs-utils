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
import groovy.util.logging.Log

/**
 * A simple tool for downsampling paired end FASTQ to achieve lower coverage.
 * 
 * @author Simon Sadedin
 */
@Log
class DownsampleFASTQ extends ToolBase {
    
    static void main(String[] args) {
        cli('DownsampleFASTQ -r <fraction of reads to preserve> -i1 <fastq r1> -i2 <fastq r2> -o1 <fastq output1> -o2 <fastq output2>', args) {
            i1 'FASTQ R1', args:1, required: true
            i2 'FASTQ R2', args:1, required: true
            o1 'Output FASTQ R1 (optional)', args:1
            o2 'Output FASTQ R2 (optional)', args:1
            r 'Fraction of reads to preserve', args:1, required: true
        } 
    }

    @Override
    public void run() {
    
        // Ensure that either both o1 and o2 options are provided or neither.
        if ((opts.o1 && !opts.o2) || (opts.o2 && !opts.o1)) {
            throw new IllegalArgumentException("Either both o1 and o2 must be provided or neither.")
        }
    
        final double rate = opts.r.toDouble()
        // If output options are provided, use them as is. Otherwise, compute defaults.
        String o1Path = opts.o1 ?: computeDefaultOutput(opts.i1, rate)
        String o2Path = opts.o2 ?: computeDefaultOutput(opts.i2, rate)
       
        OutputStream os1 = new GZIPOutputStream(new FileOutputStream(o1Path))
        OutputStream os2 = new GZIPOutputStream(new FileOutputStream(o2Path))
        os1.withWriter { w1 ->
            os2.withWriter { w2 ->
                FASTQ.filter(opts.i1, opts.i2, w1, w2) { r1, r2 ->
                    return (new Random().nextDouble() < rate)
                }            
            }
        }

        // Optional logging to show what file paths are used.
        log.info("Wrote: " + o1Path + "," + o2Path)
     }

    /**
     * Compute a default output path based on the input file and downsample rate.
     * The computed output retains the original input file name but is placed in a subdirectory 
     * named "downsampled_<percentage>" located in the same directory as the input file (or the current 
     * directory if none exists).
     */
    private String computeDefaultOutput(String inputFile, double rate) {
        int percentage = (int) (rate * 100)
        
        File inFile = new File(inputFile)
        File parentDir = inFile.parentFile ?: new File(".")
        File outDir = new File(parentDir, "downsampled_${percentage}")
        if (!outDir.exists()) {
            outDir.mkdirs()
        }
        
        // Keep the same input file name.
        return new File(outDir, inFile.name).getAbsolutePath()
    }
}
