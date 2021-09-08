package gngs.tools

import gngs.Abort
import gngs.ProgressCounter
import gngs.ToolBase
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import htsjdk.samtools.util.BlockCompressedInputStream

/**
 * Splits FASTQ in block-gzipped files which can be written live simultaneously
 * while this tool is streaming them out.
 * <p>
 * Since there is no intrinsic end point when the input FASTQs could be finished,
 * the end point is flagged by a "stop file" that must be given by the user.
 * The natural way to do that would be to put a "touch" command after the commands
 * writing the FASTQ files finishes (eg: bcl2fastq) or alternatively use some file 
 * that automatically is written at the end of the process.
 * 
 * @author Simon Sadedin
 */
@Log
class SplitLiveFASTQ extends ToolBase {
    
    static void main(String[] args) {
        cli('SplitLiveFASTQ -r1 <read1> -r2 <read2>', 'Splits FASTQ files written live from bcl2fastq to interleaved fastq while downsampling', args) {
            r1 'Read 1 BGZipped fastq', args:1, required: true
            r2 'Read 2 BGZipped fastq', args:1, required: true
            s 'Sharding specification in form n,m n=shard to output, n=total shards', args:1, required: true
            f 'file to signal input data has ended', args:1, required: true, type: File
        }
    }
    
    boolean done

    @Override
    public void run() {
        
        File stopFile = opts.f
        new Thread({
            while(true) {
                if(stopFile.exists()) {
                    log.info "Done."
                    done = true
                    return
                }
                Thread.sleep(5000)
            }
        }).start()
        
        try {
            mergeInputData()
        }
        catch(Abort abort) {
            log.info "Finished copying input files"
        }
    }

    @CompileStatic
    private mergeInputData() {
        
        List shardParts = ((String)opts['s']).tokenize(',')*.toInteger()
        int shardNumber = shardParts[0]
        int totalShards = shardParts[1]
        
        log.info "Writing shard $shardNumber out of $totalShards"
        
        List<String> r1 = new ArrayList(4)
        List<String> r2 = new ArrayList(4)

        BlockCompressedInputStream i1 = new BlockCompressedInputStream(new File((String)opts['r1']))
        BlockCompressedInputStream i2 = new BlockCompressedInputStream(new File((String)opts['r2']))

        Reader reader1 = new InputStreamReader(i1)
        Reader reader2 = new InputStreamReader(i2)

        // Read 4 lines at a time first from R1 then from R2
        ProgressCounter counter = new ProgressCounter(log:log)
        int count = 0
        while(true) {

            counter.count()
            boolean output = (count % totalShards == shardNumber)
            
            // R1
            for(int i=0; i<4; ++i) {
                String line = readPatiently(reader1)
                if(output)
                    println(line)
            }
            
            // R2
            for(int i=0; i<4; ++i) {
                String line = readPatiently(reader2)
                if(output)
                    println(line)
            }
            ++count
        }
        counter.end()
    }
    
    @CompileStatic
    String readPatiently(Reader reader) {
        while(true) {
            String line = reader.readLine()
            if(line)
                return line
                
             if(done)
                 throw new Abort()
             Thread.sleep(50)
        }
    }
}
