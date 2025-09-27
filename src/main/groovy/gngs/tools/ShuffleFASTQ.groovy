package gngs.tools

import java.util.concurrent.LinkedBlockingDeque

import gngs.FASTQ
import gngs.FASTQRead
import gngs.ToolBase
import gngs.Utils
import groovy.transform.CompileStatic

class ShuffleFASTQ extends ToolBase {

    @Override
    public void run() {
        
        String r1 = opts.r1
        String r2 = r1.replaceAll('_R1', '_R2')
        
        String o1 = r1.replaceAll('.fastq.gz$', '.shuffled.fastq.gz')
        String o2 = r1.replaceAll('.fastq.gz$', '.shuffled.fastq.gz')
        
        FASTQRead [][] buffer = new FASTQRead[(int)(opts.n)]
         writeFASTQ(buffer, r1, r2, o1, o2)
    }


    @CompileStatic
    private writeFASTQ(final FASTQRead [][] buffer, String r1, String r2, String o1, String o2) {
        Utils.withWriters([o1,o2]) {  w1, w2 ->
            FASTQ.eachPair(r1, r2) { read1, read2 ->
                // Generate random index in buffer
                int randomIndex = (int)(Math.random() * buffer.length)
                
                // If there's already a pair at this position, write it out
                if(buffer[randomIndex] != null) {
                    w1.write(buffer[randomIndex][0].toString())
                    w2.write(buffer[randomIndex][1].toString())
                }
                
                // Store the current pair at the random position
                buffer[randomIndex] = [read1, read2] as FASTQRead[]
            }
            
            // After processing all pairs, write out remaining buffered reads
            for(FASTQRead[] pair : buffer) {
                if(pair != null) {
                    w1.write(pair[0].toString())
                    w2.write(pair[1].toString()) 
                }
            }
        }
    }
    
    static void main(String[] args) {
        cli('-i1 <FASTQ>', args) {
            i1 'FASTQ R1 file to shuffle: R2 is inferred', args:1, required: true
            n 'Number of reads to include in shuffle buffer', args:1, required: true, type: Integer
        }
    }
}
