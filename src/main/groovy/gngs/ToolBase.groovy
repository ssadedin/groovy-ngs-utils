package gngs

import java.lang.invoke.MethodHandles
import java.lang.reflect.Constructor

abstract class ToolBase {
    
    OptionAccessor opts
    
    static void cli(String usage, String [] args, Closure specBuilder) {
        
        Utils.configureSimpleLogging()
        
        Cli cli = new Cli(usage: usage)
        cli.h 'Show help', longOpt: 'help' 
        
        Class originalDelegate = specBuilder.delegate
        
        def err = System.err
        
        err.println("=" * 80)
        err.println "\n" + originalDelegate.name.replaceAll('^.*\\.','') + "\n"
        err.println("=" * 80)
        err.println ""
        
        specBuilder.delegate = cli
        specBuilder()
        
        OptionAccessor opts = cli.parse(args) 
        if(!opts) {
            err.println ""
            err.println "This tool is built with Groovy NGS - the Groovy way to work with NGS data. "
            System.exit(1)
        }
            
        ToolBase tool = originalDelegate.newInstance()
        tool.opts = opts
        tool.run()
    }
    
    /**
     * Print a formatted error message and exit
     */
    void error(String msg) {
        System.err.println ""
        System.err.println "ERROR: " + msg
        System.err.println ""
        System.exit(1)
    }
    
    abstract void run() 

}
