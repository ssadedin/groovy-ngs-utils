package gngs.tools

import gngs.ProgressCounter
import gngs.SAM
import gngs.ToolBase
import groovy.util.logging.Log
import groovyx.gpars.GParsPool

@Log
class ShrinkBAM extends ToolBase {

    @Override
    public void run() {
        
        if(!opts.o.exists()) {
            throw new IllegalArgumentException("Could not access specified output directory $opts.o")
        }
        
        if(!opts.arguments()) {
            throw new IllegalArgumentException("Please specify one or more BAM files to shrink")
        }
        
        int readNameLength = opts.l?:10
        
        GParsPool.withPool(opts.n?:2) {
            
            opts.arguments().eachParallel { bam ->
                
                ProgressCounter progress = new ProgressCounter(withTime:true, log:log)
                String outPath = "$opts.o/${new File(bam).name}"
                log.info "Shrinking $bam to $outPath ..."
                new SAM(bam).filter(outPath) { r ->
                    r.baseQualities = []
                    r.readName = gngs.Hash.sha1(r.readName)[0..10]
                    return true
                }
                log.info "Indexing $outPath"
                SAM.index(new File(outPath))
                log.info "Finished $outPath"
            }
        }
    }
    
    static void main(String [] args) {
        cli('ShrinkBAM [options] <bam files...>', 'Shrinks BAM files by removing quality scores and shortening read names', args) {
            o 'Output directory', args:1, type: File, required: true
            l 'Read name length (10, use longer if read names become degenerate)', longOpt: 'read-name-length', args:1, type: Integer, required: false
            n 'Threads to use', args:1, type: Integer
        }
    }

}
