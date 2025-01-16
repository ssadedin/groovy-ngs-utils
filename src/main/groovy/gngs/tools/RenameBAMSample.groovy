package gngs.tools

import gngs.SAM
import gngs.ToolBase
import groovy.transform.CompileStatic
import groovy.util.logging.Log

/**
 * Simple utility to rename a sample in a BAM file and indexing the new file
 * 
 * (can be done by reheadering a BAM, so this is a only for convenience)
 */
@Log
class RenameBAMSample extends ToolBase {
    
    @Override
    @CompileStatic
    public void run() {
        String inputBAM = opts.arguments()[0]
        String outputSample = opts['s']
        
        File outputFile = (File)opts['o']
        
        SAM bam = new SAM(inputBAM)
        if(bam.samples.size() > 1) 
            throw new IllegalArgumentException("This tool only works with single sample BAM files, but the provided BAM file has samples: ${bam.samples.join(',')}")
        
        log.info "Renaming sample ${bam.samples[0]} in $inputBAM to $outputSample"
        bam.withWriter(sampleId: outputSample, createIndex: true, outputFile.absolutePath, true) {  w ->
            bam.eachRecord { r ->
                w.addAlignment(r)
            }
        }
        
        File indexFile = new File(outputFile.absoluteFile.parentFile.path, outputFile.name.replaceAll('.bam$', '.bai'))
        if(indexFile.exists())
            indexFile.renameTo(new File(outputFile.path + '.bai'))
        
        log.info "Wrote ${opts['o']}"
    }
    
    static void main(String [] args) {
        cli('RenameBAMSample -s <bam file>', args) {
            s 'New sample id', longOpt: 'sample', args:1, type: String, required: true
            o 'Output BAM file path', longOpt: 'output', args:1, type: File, required: false
        }
    }
}
