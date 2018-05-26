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
 * 
 * @author Simon Sadedin
 */
@Log
abstract class ToolBase {
    
    Cli parser
    
    OptionAccessor opts
    
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
        
        setProxy()
            
        ToolBase tool = originalDelegate.newInstance()
        tool.parser = cli
        tool.opts = opts
        tool.run()
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
