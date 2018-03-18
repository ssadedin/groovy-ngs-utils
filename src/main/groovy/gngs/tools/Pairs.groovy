package gngs.tools

import gngs.*
import gngs.pair.PairLocator
import gngs.pair.PairScanner
import groovy.util.logging.Log
import groovyx.gpars.actor.DefaultActor


@Log
class Pairs extends ToolBase {
    
    SAM bam

    @Override
    public void run() {
        log.info "Extracting read pairs from $opts.bam"
        bam = new SAM(opts.bam)
        
        Writer out
        if(opts.o) {
            out = new BufferedWriter(new File(opts.o).newWriter(), 2024*1024)
        }
        else {
            out = System.out.newWriter()
        }
        
        out.withWriter { 
            PairScanner scanner = new PairScanner(out, opts.n ? opts.n.toInteger():4, opts.L?getRegions():null)
            if(opts.s) {
                if(!opts.s ==~ /[0-9]*,[0-9*]/)
                    throw new IllegalArgumentException("Please provide shard number and total number of shards in form s,N to -s")
                    
                scanner.shardId = opts.s.tokenize(',')[0].toInteger()-1
                if(scanner.shardId<0)
                    throw new IllegalArgumentException("Please specify shard id > 0")
                    
                scanner.shardSize = opts.s.tokenize(',')[1].toInteger()
                if(scanner.shardId >= scanner.shardSize)
                    throw new IllegalArgumentException("Please specify shard id < number of shards ($scanner.shardSize)")
            }
            
            if(opts.namepos)
                scanner.formatter.addPosition = true
                
            scanner.scan(bam)
            
            // Debug option: dumps residual unpaired reads at end
            if(opts.du) {
                scanner.locators.each { PairLocator locator ->
                    if(!locator.buffer.isEmpty()) {
                        log.info "ERROR: Residual reads in locator: "
                        locator.buffer.each { key, value ->
                            log.info "$key: $value.r1ReferenceName:$value.r1AlignmentStart"
                        }
                    }
                }
            }
        }
    }
    
    Regions getRegions() {
        log.info "Initialising regions to scan from $opts.L"
        Regions regions
        if(opts.L.endsWith('.bed')) {
            regions = new BED(opts.L).load()
        }
        else {
            regions = new Regions()
            regions.addRegion(new Region(opts.L))
        }
        
        log.info "There are ${Utils.humanBp(regions.size())} included regions"
        
        if(opts.pad) {
            log.info "Padding regions by ${opts.pad}bp"
            regions = regions.widen(opts.pad.toInteger())
            regions = regions.reduce()
            log.info "After padding, regions span ${regions.size()}bp (${Utils.humanBp(regions.size())})"
        }
  
        return regions
    }
    
    static void main(args) {
        cli('Pairs -bam <bam> -L <regions>', args) {
            bam 'BAM file to extract read pairs from', args:1, required: true
            pad 'Amount to pad regions by (0)', args:1, required: false
            n 'Concurrency parameter (4)', args:1, required: false
            s 'Sharding factor: format <n>,<N>: output only reads belonging to shard n of N', args:1, required: false
            namepos 'Add original position to the read names', required:false
            'L' 'Regions to include reads (and mates of reads) from', longOpt: 'regions', args:1, required: false
            o 'Output file', args:1, required: false
        }
    }
}
