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
import java.util.zip.GZIPInputStream

import gngs.*
import groovy.transform.CompileStatic
import groovy.util.logging.Log

/**
 * Utility to filter 10X barcode off the R1 FASTQ
 *
 * @author Damien Zammit
 */
@Log
class TrimBarcode10X extends ToolBase {

    static void main(String [] args) {
        cli('TrimBarcode10X -r1 <FASTQ R1> [-b <Allowed Barcodes gz file>]', args) {
            r1 'FASTQ file read 1', args: 1, required: true
            b 'Allowed barcodes gz file', args: 1, required: true
        }
    }

    @Override
    public void run() {
        Writer gzOutput = new BufferedOutputStream(new GZIPOutputStream(System.out), 2048*1024).newWriter()
        try {
            run(gzOutput)
        }
        finally {
            gzOutput.close()
        }
    }

    public Map parse10XWhitelist(String file) {
        def barcodes = [:]

        if(file == "") {
            return barcodes
        }

        BufferedInputStream barcode10XWhitelistStream = new BufferedInputStream(new GZIPInputStream(new FileInputStream(file)))
        try {
            barcode10XWhitelistStream.withReader { Reader reader ->
                while(true) {
                    String allowedBarcode = reader.readLine()
                    if(allowedBarcode == null)
                        break
                    barcodes.put(allowedBarcode, true)
                }
            }
        }
        finally {
            barcode10XWhitelistStream.close()
        }
        return barcodes
    }

    @CompileStatic
    public run(Writer output) {
        int numberProcessed = 0
        Map allowedBarcodes = parse10XWhitelist((String)opts['b'])
        FASTQ.filter10Xnobarcode((String)opts['r1'], output) { FASTQRead r1Rest, FASTQRead barcode10X ->
            String barcode10Xcmp = barcode10X.bases.drop(1)
            if (allowedBarcodes.containsKey('A' + barcode10Xcmp)) {
                barcode10X.bases = barcode10X.bases + '-1,A' + barcode10Xcmp
            } else if (allowedBarcodes.containsKey('C' + barcode10Xcmp)) {
                barcode10X.bases = barcode10X.bases + '-1,C' + barcode10Xcmp
            } else if (allowedBarcodes.containsKey('G' + barcode10Xcmp)) {
                barcode10X.bases = barcode10X.bases + '-1,G' + barcode10Xcmp
            } else if (allowedBarcodes.containsKey('T' + barcode10Xcmp)) {
                barcode10X.bases = barcode10X.bases + '-1,T' + barcode10Xcmp
            } else {
	        r1Rest.bases = barcode10X.bases + r1Rest.bases
		r1Rest.quals = barcode10X.quals + r1Rest.quals
	    }
            ++numberProcessed
            return true
        }
        log.info "Wrote ${numberProcessed}"
    }
}
