package gngs.tools

import gngs.*
import groovy.transform.CompileStatic
import groovy.util.logging.Log

/**
 * Utility class for tracking state of a specific VCF in the merge process
 * 
 * @author Simon Sadedin
 */
@CompileStatic
class VCFMergeSource implements Comparable {

    Variant next
    Iterator<Variant> iter

    @Override
    public int compareTo(Object o) {
        VCFMergeSource ov = (VCFMergeSource)o
        return this.next?.xpos <=> ov.next?.xpos
    }
}

/**
 * Merges together multiple sorted VCFs that contain the same samples with non-overlapping regions. Designed
 * to facilitate scatter-gather type concurrency. Similar to Picard GatherVCFs but 
 * that program requires VCFs to be strictly ordered, whereas this program allows them to 
 * have interleaveed regions. Also similar to {@link MergePartitionedVCFs}, but this class is guaranteed
 * to use minimal memory and have single-threaded stable performance, while MergePartitionedVCFs may perform better on
 * large numbers of VCFs because it parallelises reads but at the cost of requiring memory for buffering.
 * 
 * @author Simon Sadedin
 */
@Log
class MergeInterleavedVCFs extends ToolBase {

    @Override
    public void run() {
        
        log.info "Merging ${opts.is.size()} VCFs from " + opts.is.join(', ')

        List<VCF> vcfs = opts.is.collect { new VCF(it) }
        TreeSet<VCFMergeSource> nextVariants = new TreeSet<VCFMergeSource>()
        
        vcfs.each { VCF vcf ->
            Iterator<Variant> iter = vcf.iterator()
            if(iter.hasNext())
                nextVariants.add(new VCFMergeSource(iter:iter, next: iter.next()))
        }

        Utils.writer(opts.o).withWriter { w ->
            vcfs[0].printHeader(w)
            writeVariants(nextVariants, w)
        }
        
        log.info "Wrote ${opts.o}"
        log.info "Done"
    }
    
    @CompileStatic
    void writeVariants(final TreeSet<VCFMergeSource> nextVariants, final Writer w) {
        Variant last = null
        ProgressCounter progress = new ProgressCounter(withRate:true, withTime:true, extra: { last.toString() })
        while(!nextVariants.isEmpty()) {
            VCFMergeSource output = nextVariants.pollFirst()
            last = output.next
            w.write(output.next.line)
            w.write('\n')
            if(output.iter.hasNext()) {
                output.next = output.iter.next()
                nextVariants.add(output)
            }
            progress.count()
        }
    }
    
    public static void main(String [] args) {
        cli('MergeInterleavedVCFs -i <vcf1> -i <vcf2> ... -o <vcf>', args) {
            i 'Input VCF', args: '+', required: true
            o 'Output VCF', args: 1, required:false
        }
    }
}
