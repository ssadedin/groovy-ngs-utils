/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2018 Simon Sadedin, ssadedin<at>gmail.com
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
package gngs

import java.lang.invoke.MethodHandles
import java.lang.reflect.Constructor

import groovy.util.logging.Log

/**
 * A support class that makes implementing a tool based on GNGS very easy 
 * while supporting some standard behaviors (such as logging, etc).
 * <p>
 * To implement a tool, use the following steps:
 * <li>Create a class extending <code>ToolBase</code>
 * <li>Implement the <code>static void main(String [] args)</code> method
 * <li>In the main method, add a call to the {@link #cli} support method,
 *     which will run the context of a {@link gngs.Cli} class to allow you
 *     to configure the command line behavior
 * <li>
 * <p>
 * A simple example:
 * 
 * <pre>
 * @Log
 * class FooTool extends ToolBase {
 *     void run() {
 *         // ... do stuff with optsion via provided <code>opts</opts> variable ...
 *         if(opts.i) {
 *             // ...
 *         }
 *     }
 * 
 *     static void main(String [] args) {
 *         cli('Does some foo stuff', args) {
 *            i 'Input to the foo', longOpt: 'input', args:1
 *         }
 *     }
 * }
 * </pre>
 * Creating tests for tools based on ToolBase encounters some challenges due to the static methods
 * and the way the command line options are translated directly to the {@link opts} variable. To help, a 
 * dedicated {@link #test} method is provided that allows invocation of the tool using command line arguments,
 * which correctly sets up the {@link opts} variable as appropriate.
 * <p>
 * Example:
 * <pre>
 *     MultiCov mc = new MultiCov()
 *     mc.test(['-bed', 'src/test/data/multicov_test.bed', 'src/test/data/small.test.bam']) {
 *         mc.run()
 *     }
 * </pre>
 *
 * 
 * @author Simon Sadedin
 */
@Log
abstract class ToolBase {
    
    Cli parser
    
    /**
     * The options parsed for the tool
     */
    CliOptions opts
    
    static ThreadLocal<ToolBase> test = new ThreadLocal<ToolBase>()
    
    static String footer = "This tool is built with Groovy NGS - the Groovy way to work with NGS data.\n" 
    
    /**
     * A utility method to facilitate testing of tools that extend this class.
     * <p>
     * To use this method, create an instance of the tools, then invoke the test method with the 
     * command line arguments to pass, and then execute actual functions within the closure callback.
     * 
     * @param args  the command line args
     * @param c
     */
    void test(List<String> args, Closure c) {
        test.set(this)
        this.main(args as String[])
        c()
    }
    
    static void cli(String usage, String [] args, Closure specBuilder) {
        cli(usage, null, args, specBuilder)
    }
    
    /**
     * Create an instance of the enclosing class and call its {@link #run} method after 
     * parsing options with a {@link gngs.Cli} instance configured by 
     * the provided <code>specBiulder</code> closure.
     * 
     * @param usage         a string to present to the user explaining usage of the tool
     * @param args          raw args from command line
     * @param specBuilder   a closure to configure a {@link groovy.util.CliBuilder} via a 
     *                      {@link gngs.Cli} instance
     */
    static void cli(String usage, String header, String [] args, Closure specBuilder) {
        
        Utils.configureSimpleLogging()
        
        if(header) {
            // lets us format it however we want in our groovy code and shows up reasonably
            // wrapped from commons-cli
            header = '\n' + header.stripIndent().readLines().grep { it }.join(' ') + '\n\nOptions:\n\n'
        }
        
        Cli cli = new Cli(usage: usage, header: header)
        cli.h 'Show help', longOpt: 'help' 
        
        Class originalDelegate = specBuilder.delegate
        
        def err = System.err
        specBuilder.delegate = cli
        specBuilder()
        
        if('-h' in args || '--help' in args) {
            printHeader(originalDelegate, null)
            printHelpAndExit(cli, err)
        }
        
        OptionAccessor rawOpts = cli.parse(args) 
        if(!rawOpts) {
            printHeader(originalDelegate, null)
            err.println ""
            err.println footer
            System.exit(1)
        }
        CliOptions opts = new CliOptions(opts:rawOpts)
        setProxy()
            
        ToolBase tool = test.get()?: originalDelegate.newInstance()
        tool.parser = cli
        tool.opts = opts
        
        printHeader(originalDelegate, tool)
        
        if(test.get()) 
            return
            
        try {
            tool.run()
        }
        catch(IllegalArgumentException e) {
            System.err.println "One or more arguments were invalid or missing: " + e.message + "\n"
            cli.usage()
            err.println "\n" + footer 
            System.exit(1)
        }
    }
    
    private static printHeader(Class originalDelegate, ToolBase tool) {
        def err = System.err 
        err.println("=" * 80)
        err.println("")
        if(tool)
            tool.printTitle()
        else
            System.err.println originalDelegate.name.replaceAll('^.*\\.','') 
        err.println("")
        err.println("=" * 80)
        err.println ""
    }
    
    void printTitle() {
        System.err.println this.class.name.replaceAll('^.*\\.','') 
    }

    private static printHelpAndExit(Cli cli, PrintStream err) {
        cli.usage()
        err.println ""
        err.println footer
        System.exit(0)
    }
    
    /**
     * Check if the http_proxy variable is set, and if no proxy set, initialize the system properties
     * needed from there.
     */
    static void setProxy() {
        String envProxy = System.getenv('http_proxy')
        if(envProxy && !System.properties.containsKey('http.proxyHost')) {
            // Parse the proxy out
            
            String host = null
            String port = null
            List<String> parts = envProxy.tokenize(':')
            host = parts[1].replaceAll('^//','')
            if(parts.size()>2) { // Port set
                port = parts[2]
            }
            log.info "Auto detected proxy host=$host, proxy port=$port"
            System.properties['http.proxyHost'] = host
            
            if(port != null)
                System.properties['http.proxyPort'] = port
        }
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
