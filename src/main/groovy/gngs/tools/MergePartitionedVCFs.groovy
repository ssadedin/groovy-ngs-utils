package gngs.tools

import java.util.concurrent.atomic.AtomicInteger

import gngs.*
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import groovyx.gpars.GParsPool




@Log
class MergePartitionedVCFs extends ToolBase {
    
    VCFMergerActor merger
    
    List<String> contigs = null

    @Override
    public void run() {
        
        VCF headerVCF = new VCF(opts.i)
        this.contigs = headerVCF.contigs as String[]
        
        log.info "Merging ${opts.is.join(',')}"
        Utils.writer(opts.o).withWriter { w ->
            merger = new VCFMergerActor(headerVCF, opts.is.size(), w)
            merger.progress = new ProgressCounter(withRate: true, withTime:true, log:log)
            merger.start()
            
            GParsPool.withPool(opts.is.size()) {
                opts.is.eachWithIndexParallel { String vcfPath, int index ->
                    readVCF(vcfPath, index)
                }
            }
            
            log.info "Finished writing: sending stop message"
            merger.sendStop()
            merger.join()
            log.info "Stopped."
        }
        log.info "Wrote $opts.o"
    }
    
    @CompileStatic
    void readVCF(String vcfPath, int index) {
        
        AtomicInteger counter = new AtomicInteger(0)
        
        VCF.parse(vcfPath, progress:false) { Variant v -> 
            final int contigIndex = contigs.indexOf(v.chr)
            assert contigIndex >= 0
            MergeVariant m = new MergeVariant(contigIndex, v, index)
            merger.sendLimited(new AcknowledgeableMessage(m, counter))
            return false
        }
    }
    
    public static void main(String [] args) {
        cli('MergePartitionedVCFs -i <vcf1> -i <vcf2> ... -o <vcf>', args) {
            i 'Input VCF', args: Cli.UNLIMITED, required: true
            o 'Output VCF', args: 1, required:false
        }
    }
}

@CompileStatic
class MergeVariant implements Comparable<MergeVariant> {
    int contig
    Variant  variant
    int sourceIndex
    
    public MergeVariant(int contig, Variant variant, int sourceIndex) {
        super();
        this.contig = contig;
        this.variant = variant;
        this.sourceIndex = sourceIndex
    }


    @Override
    public int compareTo(MergeVariant v2) {
        int result = contig - v2.contig;
        if(result)
            return result
        return variant.pos - v2.variant.pos
    } 
}

@Log
class VCFMergerActor extends RegulatingActor<MergeVariant> {
    
    VCF header
    
    Writer out
    
    final int numVCFs
    
    PriorityQueue<MergeVariant> queue = new PriorityQueue(100000)
    
    Integer [] counts 
    
    VCFMergerActor(VCF headerVcf, int numVCFs, Writer out) {
        super(100, 1000)
        this.numVCFs = numVCFs
        this.counts = new Integer[numVCFs]
        for(int i=0; i<numVCFs; ++i) {
            this.counts[i] = 0i
        }
        
        this.out = out
        out.println(headerVcf.headerLines.join('\n'))
    }

    @CompileStatic
    @Override
    public void process(MergeVariant v) {
        int prevContig = -1
        
        queue.add(v)
        
        if(this.counts[v.sourceIndex] == 0) {
            
            ++this.counts[v.sourceIndex]
            while(this.counts.min() > 0) {
                MergeVariant top = queue.poll()
                assert top.contig >= prevContig
                out.write(top.variant.line)
                out.write('\n')
                prevContig = top.contig
                --this.counts[top.sourceIndex]
            }
        }
        else {
            ++this.counts[v.sourceIndex]
        }
    }

    @Override
    public void onEnd() {
        log.info "Flushing ${queue.size()} residual entries ..."
        while(!queue.isEmpty()){
            out.println(queue.poll().variant.line)
        }
    }
    
    @CompileStatic
    @Override
    String toString() {
        "VCFMergerActor(${numVCFs} vcfs)"
    }
}