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
 * Utility to split paired end FASTQ and write it out in interleaved 10X format
 *
 * @author Damien Zammit
 */
@Log
class SplitFASTQ10X extends ToolBase {

    static void main(String [] args) {
        cli('SplitFASTQ10X -s <Sharding Specification> -r1 <FASTQ R1> -r2 <FASTQ R2> -i1 <FASTQ I1> [-b <Allowed Barcodes gz file>]', args) {
            s 'Sharding specification in form n,m n=shard to output, n=total shards', args:1, required: true
            r1 'FASTQ file read 1', args: 1, required: true
            r2 'FASTQ file read 2', args: 1, required: true
            i1 'FASTQ file index 1', args: 1, required: true
            b 'Allowed barcodes gz file', args: 1, required: false
        }
    }

    @Override
    public void run() {
        List shardParts = opts.s.tokenize(',')
        int shardId = shardParts[0].toInteger()
        int numberOfShards = shardParts[1].toInteger()

        log.info "Splitting ${opts.r1},${opts.r2} ${numberOfShards} ways and outputting shard ${shardId}"

        Writer gzOutput = new BufferedOutputStream(new GZIPOutputStream(System.out), 2048*1024).newWriter()
        try {
            run(shardId, numberOfShards, gzOutput)
        }
        finally {
            gzOutput.close()
        }
    }

    public String[] parse10XWhitelist(String file) {
        def barcodes = []
        barcodes << ""

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
                    barcodes << allowedBarcode
                }
            }
        }
        finally {
            barcode10XWhitelistStream.close()
        }
        return barcodes
    }

    @CompileStatic
    public run(int shardId, int shards, Writer output) {
        int numberWritten = 0
        int numberProcessed = 0
        String[] allowedBarcodes = parse10XWhitelist((String)opts['b'])
        FASTQ.filter10X((String)opts['r1'], (String)opts['r2'], (String)opts['i1'], output) { FASTQRead r1Rest, FASTQRead r2, FASTQRead barcode10X, FASTQRead i1 ->
            int hash = r2.name.hashCode()
            int readShardId = Math.abs(hash % shards)
            ++numberProcessed
            if(readShardId == shardId) {
                if (allowedBarcodes == [""] || barcode10X.bases in allowedBarcodes) {
                    ++numberWritten
                    return true
                }
            }
            else {
               return false
            }
        }
        log.info "Wrote ${numberWritten} / ${numberProcessed} (${Utils.perc(numberWritten/numberProcessed)}) for ${shards} shards"
    }
}
