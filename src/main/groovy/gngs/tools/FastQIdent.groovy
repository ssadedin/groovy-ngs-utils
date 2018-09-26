package gngs.tools

import java.nio.charset.Charset

import com.google.common.hash.BloomFilter
import com.google.common.hash.Funnels

import gngs.*
import groovy.transform.CompileStatic
import groovy.util.logging.Log

/**
 * Find sequences not present in a raw fastq file
 * 
 * @author Simon Sadedin
 */
@Log
class FastQIdent extends ToolBase {
    
    BloomFilter<String> filter = BloomFilter.create(Funnels.stringFunnel(Charset.forName('UTF-16')), 1000000000, 0.001d )

    static void main(String [] args) {
        cli('-f1 <fastq1> -f2 <fastq1>', args) {
            f1 'First FASTQ',args:1,required: true
            f2 'Second FASTQ', args:1,required: true
        }
    }
    
    @CompileStatic
    private static String getComparableContent(FASTQRead r) {
        r.name.tokenize(' ')[0].tokenize('/')[0] + r.bases        
    }

    @Override
    @CompileStatic
    public void run() {
        
        String f1 = opts['f1']
        String f2 = opts['f2']
        
        log.info "Comparing $f1 to $f2"
        
        ProgressCounter f1Progress = new ProgressCounter(withTime: true, withRate:true)
        FASTQ.eachRead(f1) { FASTQRead r1 ->
            filter.put(getComparableContent(r1))
            f1Progress.count()
        }
        
        log.info "Now identifying missing reads ..."
        int missingReads = 0
        ProgressCounter f2Progress = new ProgressCounter(withTime: true, withRate:true, extra: {
            "Missing reads: $missingReads"
        })
        
        System.out.withWriter { w -> 
            FASTQ.eachRead(f2) { FASTQRead r1 ->
                if(!filter.mightContain(getComparableContent(r1))) {
                    System.err.println "Read $r1.name not found or not identical"
                    r1.write(w)
                }
                f2Progress.count()
            }
        }
        
        log.info "Done."
        log.info "Identified ${missingReads} reads that were non-identical in file2, compared to file 1"
    }
}
