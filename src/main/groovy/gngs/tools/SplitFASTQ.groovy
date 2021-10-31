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

import java.io.Writer
import java.nio.charset.Charset
import java.util.concurrent.atomic.AtomicInteger
import java.util.zip.GZIPOutputStream

import com.google.common.hash.BloomFilter
import com.google.common.hash.Funnels

import gngs.*
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import groovyx.gpars.actor.Actor


@Log
class FASTQDedupeActor extends RegulatingActor {

    static BloomFilter filter = BloomFilter.create(Funnels.stringFunnel(Charset.forName('UTF-16')), 800000000, 0.01)
    
    static long duplicates = 0
    static long total = 0
    
    FASTQDedupeActor(RegulatingActor downstream) {
        super(downstream, 5000,10000)
    }
    
    @Override
    @CompileStatic
    public void process(Object message) {
        
        FASTQRead [] reads = (FASTQRead[]) message
        FASTQRead r1 = reads[0]
        FASTQRead r2 = reads[1]
        
        String forward = r1.bases + r2.bases
        String reverse = FASTA.reverseComplement(forward)
        
        if(filter.mightContain(forward) || filter.mightContain(reverse)) {
            ++duplicates
        }
        else {
            filter.put(forward)
            this.sendDownstream(message)
        }        
        ++total
    }
}

@Log
@CompileStatic
class FASTQWriterActor extends RegulatingActor {
    
    Writer output
    
    int numberWritten = 0
    int numberProcessed = 0
    
    int shardId
    int shards

    public FASTQWriterActor(Writer output, int shardId, int shards) {
        super(5000, 20000);
        this.output = output;
        this.shardId = shardId
        this.shards = shards
    }
    
    @Override
    public void process(Object message) {
        FASTQRead [] reads = (FASTQRead[]) message
        
        FASTQRead r1 = reads[0]
        FASTQRead r2 = reads[0]
        
        int hash = r1.name.hashCode()
        
        int readShardId = Math.abs(hash % shards)
        ++numberProcessed
        if(readShardId == shardId) {
            ++numberWritten
            r1.write(this.output)
            r2.write(this.output)
        }
        println "wrote $r1.name"
    }
}

/**
 * A flexible utility to "reshape" paired end fastq data for processing by merging and splitting 
 * it and writing it out in interleaved format so that there is no need to store intermediate 
 * files.
 * 
 * @author Simon Sadedin
 */
@Log
class SplitFASTQ extends ToolBase {
    
    static void main(String [] args) {
        cli('SplitFASTQ -s <Sharding Specification> -r1 <FASTQ R1> -r2 <FASTQ R2>', args) {
            s 'Sharding specification in form n,m n=shard to output, n=total shards', args:1, required: false
            r1 'FASTQ file read 1, multiple files joined with comma', args: 1, required: true
            r2 'FASTQ file read 2, multiple files joined with comma', args: 1, required: false
            o1 'Output file read 1', args: 1, required: false
            o2 'Output file read 2', args: 1, required: false
            f 'Downsampling factor (1.0)', args: 1, required: false
            dedupe 'Only emit reads that are not duplicates of reads already seen'
            n 'Number of threads to use for deduping (4)', args:1, required: false
        }
    }
    
    boolean downsample
    
    double downsampleFactor = 1.0d
    
    Random downsampleRandomNumberGenerator

    @Override
    public void run() {
        List shardParts = opts.s ? opts.s.tokenize(',') : '0,1'
        int shardId = shardParts[0].toInteger()
        int numberOfShards = shardParts[1].toInteger()
        
        if(opts.f) {
            this.downsample = true
            this.downsampleFactor = opts.f.toDouble()
            this.downsampleRandomNumberGenerator = new Random()
        }
        
        log.info "Splitting ${opts.r1},${opts.r2} ${numberOfShards} ways and outputting shard ${shardId}"
        log.info "Deduplication enabled: $opts.dedupe"
        
        Writer output1 
        Writer output2 
        if(opts.o1) {
            output1 = Utils.outputWriter(opts.o1)
        }
        else {
            output1 = System.out.newWriter()
        }
        
        if(opts.o2) {
            output2 = Utils.outputWriter(opts.o2)
        }
        else {
            output2 = output1
        }

        try {
            if(opts.dedupe) {
                runDeduped(shardId, numberOfShards, (opts.n?:4).toInteger(),output1)
            }
            else {
                runNoDedupe(shardId, numberOfShards, output1, output2)
            }
        }
        finally {
            output1.close()
            if(output1 != output2)
                output2.close()
        }
    }
    
    @CompileStatic
    public void runNoDedupe(int shardId, int shards, Writer output1, Writer output2) {
        int numberWritten = 0
        int numberProcessed = 0

        List<String> r1Paths = (List<String>)((String)opts['r1']).tokenize(',')
        List<String> r2Paths = opts['r2'] ? (List<String>)((String)opts['r2']).tokenize(',') : ["false"] * r1Paths.size()
         
        [r1Paths, r2Paths].transpose().each { String r1Path, String r2Path ->
            
            log.info "Processing FASTQ pair $r1Path, $r2Path"
        
             if(r2Path != "false") {
                FASTQ.filter(r1Path, r2Path, output1, output2) { FASTQRead r1, FASTQRead r2 ->
                    int hash = r1.name.hashCode()
                    int readShardId = Math.abs(hash % shards)
                    ++numberProcessed 
                    if(this.downsample) {
                        if(this.downsampleRandomNumberGenerator.nextDouble()>this.downsampleFactor)
                            return
                    }
                    if(readShardId == shardId) {
                        ++numberWritten
                       return true 
                    }
                    else {
                       return false 
                    }
                }
            }
            else {
                FASTQ.eachRead(r1Path) { FASTQRead r1 ->
                    int hash = r1.name.hashCode()
                    int readShardId = Math.abs(hash % shards)
                    ++numberProcessed 
                    if(this.downsample) {
                        if(this.downsampleRandomNumberGenerator.nextDouble()>this.downsampleFactor)
                            return
                    } 
                    if(readShardId == shardId) {
                        ++numberWritten
                        r1.write(output1)
                    }
                }
            }
        }
        log.info "Wrote ${numberWritten} / ${numberProcessed} (${Utils.perc(numberWritten/numberProcessed)}) for ${shards} shards"
    }
    
    @CompileStatic
    public runDeduped(int shardId, int shards, int numDedupers, Writer output) {
        
        if(this.downsample)
            throw new IllegalStateException('Downsampling is not supported when deduping is enabled')
        
        FASTQWriterActor out = new FASTQWriterActor(output, shardId, shards)
        List<FASTQDedupeActor> dedupers = (1..numDedupers).collect { new FASTQDedupeActor(out) }
        
        List<Actor> actors = (List<Actor>)([out, *dedupers])
        actors*.start()
        
        AtomicInteger pending = new AtomicInteger(0)
        int i = 0
        List<String> r1Paths = (List<String>)((String)opts['r1']).tokenize(',')
        List<String> r2Paths = opts['r2'] ? (List<String>)((String)opts['r2']).tokenize(',') : ["false"] * r1Paths.size()
         
        [r1Paths, r2Paths].transpose().each { String r1Path, String r2Path ->
         
            FASTQ.eachPair(r1Path, r2Path) { FASTQRead r1, FASTQRead r2 ->
                dedupers[++i % numDedupers].sendTo([r1,r2])
            }
        }
        
        log.info "Stopping dedupers ...."
        dedupers*.sendStop()
        dedupers*.join()
        log.info "Stopping output writer ...."
        out << "stop"
        out.join()
        log.info "Actors stopped."
        
        log.info "Deduplication rate: " + FASTQDedupeActor.duplicates + " / " + FASTQDedupeActor.total + " (${Utils.perc(FASTQDedupeActor.duplicates /  FASTQDedupeActor.total)})"
    }
}
