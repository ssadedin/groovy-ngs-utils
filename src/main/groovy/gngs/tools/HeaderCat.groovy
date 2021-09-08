package gngs.tools

import gngs.ToolBase
import gngs.Utils
import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.util.logging.Log

/**
 * A utility to concatenate files that may have headers which 
 * should only be included once at the top.
 * 
 * @author Simon Sadedin
 */
@Log
@CompileStatic
class HeaderCat extends ToolBase {

    @Override
    public void run() {
        List<String> vcfPaths = opts.arguments()
        Utils.writer((String)opts['o']).withWriter { Writer w ->
            List<Reader> files = (List<Reader>)vcfPaths.collect {
                Utils.reader(it)
            }
            
            String commentChar = opts['c']?:'#'
            boolean stripFirstLineOnly = (commentChar == 'NONE')
            
            int fileIndex = 0
            String vcfPath = vcfPaths[fileIndex]
            
            log.info "Processing file $vcfPath (${fileIndex+1}/${vcfPaths.size()})"

            // For the first file cat the whole file
            Reader r0 = (Reader)files[0]
            w << r0
            
            for(Reader r in files[1..-1]) {
                vcfPath = vcfPaths[++fileIndex]
                log.info "Processing file $vcfPath"
                String line = r.readLine()
                
                if(stripFirstLineOnly)
                    line = r.readLine()
                else
                    while(line?.startsWith('#'))
                        line = r.readLine()

                if(line) {
                    w << line 
                    w.write('\n')
                }

                w << r
            }
            log.info "Done"
        }
    }
    
    @CompileDynamic
    static void main(String[] args) {
       cli('HeaderCat -o <output> <file1> [<file2>] ...', 
           'Smart concatenation of files of mixed compression and possibly containing headers',
            args)  {
           c 'Comment char, or NONE if first line should be arbitrarily removed (#)', args: 1, required: false
           o 'The output file', required: true, longOpt: 'output', args:1
       }
    }

}
