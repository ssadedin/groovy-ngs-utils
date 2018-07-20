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

import gngs.*
import groovy.transform.CompileStatic
import groovy.util.logging.Log

/**
 * Utility to split paired end FASTQ and write it out in interleaved format
 * 
 * @author Simon Sadedin
 */
@Log
class SplitFASTQ extends ToolBase {

    static void main(String [] args) {
        cli('SplitFASTQ -s <Sharding Specification> -r1 <FASTQ R1> -r2 <FASTQ R2>', args) {
            s 'Sharding specification in form n,m n=shard to output, n=total shards', args:1, required: true
            r1 'FASTQ file read 1', args: 1, required: true
            r2 'FASTQ file read 2', args: 1, required: true
        }
    }

    @Override
    public void run() {
        List shardParts = opts.s.tokenize(',')
        int shardId = shardParts[0].toInteger()
        int numberOfShards = shardParts[1].toInteger()
        
        log.info "Splitting ${opts.r1},${opts.r2} ${numberOfShards} ways and outputting shard ${shardId}"
        
        Writer output = System.out.newWriter()
        try {
            run(shardId, numberOfShards, output)
        }
        finally {
            output.close()
        }
    }
    
    @CompileStatic
    public run(int shardId, int shards, Writer output) {
        int numberWritten = 0
        int numberProcessed = 0
        FASTQ.filter((String)opts['r1'], (String)opts['r2'], output) { FASTQRead r1, FASTQRead r2 ->
            int hash = r1.name.hashCode()
            int readShardId = Math.abs(hash % shards)
            ++numberProcessed 
            if(readShardId == shardId) {
                ++numberWritten
               return true 
            }
            else {
               return false 
            }
        }
        log.info "Wrote ${numberWritten} / ${numberProcessed} (${Utils.perc(numberWritten/numberProcessed)}) for ${shards} shards"
    }
}
