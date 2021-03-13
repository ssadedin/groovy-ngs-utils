package gngs.tools

import gngs.ToolBase
import gngs.Utils
import groovy.util.logging.Log

@Log
class HeaderCat extends ToolBase {

    @Override
    public void run() {
        List<String> vcfPaths = opts.arguments()
        Utils.writer(opts.o).withWriter { Writer w ->
            def files = vcfPaths.collect {
                Utils.reader(it)
            }
            
            int fileIndex = 0
            String vcfPath = vcfPaths[fileIndex]
            
            log.info "Processing file $vcfPath"

            // For the first file cat the whole file
            Reader r0 = files[0]
            w << r0
            
            for(Reader r in files[1..-1]) {
                vcfPath = vcfPaths[++fileIndex]
                log.info "Processing file $vcfPath"
                String line = r.readLine()
                while(line?.startsWith('#'))
                    line = r.readLine()
                if(line)
                    w << line
                w << r
            }
            log.info "Done"
        }
    }
    
    static void main(String[] args) {
       cli('HeaderCat -o <output> <file1> [<file2>] ...', 
           'Smart concatenation of files of mixed compression and possibly containing headers',
            args)  {
           o 'The output file', required: true, longOpt: 'output', args:1
       }
    }

}
