package gngs.tools

import gngs.SAM
import gngs.ToolBase
import groovy.transform.CompileDynamic
import groovy.util.logging.Log
import htsjdk.samtools.SAMRecord

/**
 * Tool to count the frequency of kmers at the start
 * of each read in a BAM file. This is intended to inform
 * modeling of potential shearing bias.
 * 
 * @author simon.sadedin
 */
@groovy.transform.CompileStatic
@Log
class ShearingKmerCounter extends ToolBase {
    
    static int [] translation = new int[100]
    
    static int [] reverse_translation = new int[4]
    
    final private static int KMER_SIZE = 5
    
    static {
        // The translation stores the int encodings of 'A','C','T','G'
        // in a sparse array so that we avoid needing to do any type of has lookup
        translation[65] = 0
        translation[67] = 1
        translation[71] = 2
        translation[84] = 3
        
        reverse_translation[0] = 'A' as char
        reverse_translation[1] = 'C' as char
        reverse_translation[2] = 'G' as char
        reverse_translation[3] = 'T' as char
    }
    
    /**
     * Convert a kmer encoded as an integer back to string form
     */
    String decodeKmer(int kmer_int) {
        def result = []
        for(int i=0; i<KMER_SIZE; ++i) {
            result.add(reverse_translation[kmer_int>>(i*2) & 0x3])
        }
        String kmer = result*.asType(char).join('')
        return kmer
    }
    
    static int [] countMotifs(SAM bam) {
        return (int[])bam.withIterator {
            countMotifs(it)
        }
    }
    
    static int [] countMotifs(htsjdk.samtools.SAMRecordIterator iter) {
        final int [] counts = new int[1 <<(KMER_SIZE*2)]
        while(iter.hasNext()) {
            SAMRecord read = iter.next()
            int val = 0
            final byte [] bases = read.readBases
            // If read is to the left  ---> then take first bases
            if(read.alignmentStart <= read.mateAlignmentStart) {
                for(int i=0; i<KMER_SIZE; ++i) {
                    val |= (translation[bases[i]] << (i*2))
                }
            }
            else {
                // this read is to the "right" of the other read
                // so we want the bases at the end
                // ---r-->  <---r---
                final int end = bases.size() - 1
                for(int i=0; i<KMER_SIZE; ++i) {
                    val |= (translation[bases[end-i]] << (i*2))
                }
            }
            ++counts[val]
        }
        return counts
    }
   
    static int computeKmer(SAMRecord read) {
        int val = 0
        final int KMER_SIZE=5
        final byte [] bases = read.readBases
        for(int i=0; i<KMER_SIZE; ++i) {
            val |= (translation[bases[i]] << (i*2))
        }
        return val
    }

    @Override
    @CompileDynamic
    public void run() {
        log.info "Processing $opts.i"
        SAM bam = new SAM(opts.i)
        int [] result = countMotifs(bam)
        opts.o.withWriter { Writer w ->
            w.write(
                ['sample', *(0..1023).collect { decodeKmer(it)}].join('\t')
            )
            w.write('\n')
            w.write(bam.samples[0])
            w.write('\t')
            w.write(result.join('\t'))
            w.write('\n')
        }
        log.info "Finished, wrote $opts.o"
    }
    
    @CompileDynamic
    static void main(String[] args) {
        cli('ShearingKmerCounter -i <bam file> -o <output file>', 'Counts kmers at the beginning of reads', args) {
            i 'BAM file to use', longOpt: 'input', args:1, required: true, type: File
            o 'Output file', longOpt: 'output', args:1, required: true, type: File
        }
    }
}

