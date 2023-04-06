package gngs
/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2013 Simon Sadedin, ssadedin<at>gmail.com
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

import java.util.Iterator;
import java.util.zip.GZIPInputStream

import groovy.transform.CompileStatic;

/**
 * Support for parsing, filtering, and manipuulating BED files.
 * <p>
 * This class provides and implementation for {@link Regions} that is backed
 * by a file in BED format.
 * <p>
 * Some methods hit the file directly, while others rely on reading it 
 * into memory. In cases where the data set will not fit into memory, you should
 * check how the method is going to operate and choose accordingly.
 * <p>
 * Note that the constructor does <i>not</i> load the contents of the BED
 * file. To actually load the BED file you need to call {@link #load()}.
 * <p>
 * The fourth column of the BED file is referred to as the <code>extra</code>
 * attribute. This is not loaded by default. You can request that it be loaded
 * by adding the option <code>withExtra</code> to either the constructor or
 * the {@link #load()} method.
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class BED extends Regions {
    
    /**
     * Some functions require loading of the BED file into memory. Others can use 
     * it optionally. So keep a flag to know if the BED was loaded, and if it was,
     * use memory for all operations to be faster.
     */
    boolean isLoaded = false
    
    /**
     * Use 'extra' field to allow additional data to be attached to ranges in the bed file
     */
    boolean withExtra = false
    
    /**
     * The stream from which bed line will be read
     */
    InputStream bedFileStream = null
    
    /**
     * File from which bed is read
     */
    File bedFile
    
    /**
     * Columns for reading the position data - these can be overridden to allow
     * loading of any "bed-like" file
     */
    int chrColumn = 0
    int startColumn = 1
    int endColumn = 2
    
    boolean readColumnNames=false
    
    /**
     * If true, column names will be set as properties on the regions
     */
    boolean withProperties=false
    
    List<String> columnNames = null
    
    /**
     * Empty bed file
     */
    BED(Map attributes=[:]) {
        if(attributes && attributes.withExtra != null)
            this.withExtra = attributes.withExtra
    }
    
    BED(Map attributes=[:], String fileName) {
        this(attributes, new File(fileName))
    }
    
    BED(Map attributes=[:], File file, Closure c = null) {
        this(attributes)
        this.bedFile = file
        this.bedFileStream = new FileInputStream(file)
    }
    
    BED(Map attributes=[:], InputStream inStream, Closure c = null) {
        this(attributes)
        this.bedFileStream = inStream
    }
    
    /**
     * Create a BED from a list of regions
     */
    BED(Map attributes=[:], Iterable<Region> regions) {
        this(attributes)
        for(Region r in regions) {
            addRegion(r.chr, r.from, r.to)
        }
    }
   
    /**
     * Iterate over each line in the BED file and invoke the specified
     * closure with bed file information for the line.
     * <p>
     * Options can be provided as a Map argument, including:
     * 
     * <li>unique (whether identical ranges should be included multiple times or only once)
     */
    @CompileStatic
    void eachRange(Map options, Closure c) {
        
        if(bedFileStream == null || isLoaded) { // Used to be && isLoaded - why?
            super.eachRange(options,c)
            return
        }
        
        int count = 0
        def includeInfo = c.getMaximumNumberOfParameters()>3;
        
        try {
            if(bedFileStream == null || bedFileStream.available()<0) {
                if(bedFile != null) {
                  if(bedFile.name.endsWith(".gz")) {
                      bedFileStream = new GZIPInputStream(new FileInputStream(bedFile))
                  }
                  else
                      bedFileStream = new FileInputStream(bedFile)
                }
                else
                try { bedFileStream.reset() } catch(Exception e) {}
            }
        }
        catch(IOException ex) {
            try { bedFileStream.reset() } catch(Exception e) {}
            if(bedFile != null) 
                  bedFileStream = new FileInputStream(bedFile)
         }
        
        if(bedFileStream == null)
            return
        
        HashSet processed = new HashSet()
        boolean unique = options.unique?true:false
        
        c.delegate = new Object() { void abort() { throw new Abort() } }
        
        ProgressCounter progress
        if(options.progress == true)
            progress = new ProgressCounter(lineInterval:10, timeInterval:10000)
        else
        if(options.progress instanceof ProgressCounter) {
            progress = (ProgressCounter)options.progress // User can set a progress object
            c.delegate = progress
        }
            
        try {
          bedFileStream.eachLine { String line ->
              
              ++count
              if(line.startsWith('#'))
                  return
              if(line.startsWith("track"))    
                  return
              if(line.startsWith("browser"))    
                  return
                  
              List fields = line.tokenize('\t')
              
              if(count == 1 && readColumnNames) {
                  columnNames = fields 
                  return
              }
              
              if(unique) {
                  String key = fields[chrColumn]+':'+fields[startColumn]+':'+fields[endColumn]   
                  if(processed.contains(key)) {
                      if(progress)
                          progress.count()
                      return
                  }
                  processed.add(key)
              }
              
              String chr = fields[chrColumn]
              int start = -1, end = -1
              try {
                start = Integer.parseInt(fields[startColumn])
                end = Integer.parseInt(fields[endColumn])
              }
              catch(Exception e) {
                  System.err.println "Unable to parse line $count : \n$line\n\nException reported: $e"   
                  throw e
              }
              
             try {
                if(includeInfo==true)
                    c.call(chr,start,end,fields.size()>3?fields[3]:null)
                else {
                    c.call(chr,start,end)
                }
                
                if(progress)
                    progress.count()
              }
              catch(Exception e) {
                  System.err.println "Failure while handling line $count : \n$line\n\nException reported: $e"   
                  throw e
              }
           }
        }
        catch(Abort e) {
            System.err.println "Iteration aborted at line $count"
        }
        finally {
            if(progress)
                progress.end()
        }
    }
    
    void eachLoadedRange(Closure c) {
        super.eachRange([:],c)
    }
    
//    @CompileStatic
    void eachLoadedRange(Map options, Closure c) {
       super.eachRange(options,c) 
    }
   
    /**
     * Call closure c for each range that overlaps the given position
     */
    @CompileStatic
    void eachOverlap(String chr, int pos, Closure c) {
        if(!isLoaded) {
          eachRange([:]) { String rangeChr, int rangeStart, int rangeEnd ->
              if(rangeChr == chr && pos >= rangeStart && pos << rangeEnd) {
                  c(rangeChr, rangeStart, rangeEnd)
              }
          }
        }
        else {
			super.eachOverlap(chr, pos, c)
        }
    }
   
    /**
     * Iterate each position in a BED file and call the given closure with it
     * <p>
     * NOTE: this only works on an actual BED file, not in-memory BED file.
     * 
     * @TODO: merge overlapping regions
     */
    @CompileStatic 
    void eachPosition(Closure c) {
        long totalSize = super.size()
        System.err.println("Scanning $totalSize positions")
        ProgressCounter counter = new ProgressCounter()
          bedFileStream.eachLine { String line ->
              if(line.startsWith('#'))
                  return
                  
              String [] fields = line.split('\t')
              String chr = fields[0]
              int start = Integer.parseInt(fields[1])
              int end = Integer.parseInt(fields[2])
              for(int i=start; i<end; ++i) {
                  c(chr, i)
                  counter.count { int count ->
                      System.err.println "Processed $count / $totalSize (${count/(float)totalSize}%)"
                  }
              }
              
        }
    }
    
    /**
     * Load the data for this BED file into memory.
     * <p>
     * If withExtra is false, only the range data will be loaded. This is very compact and fast.
     * Otherwise, the optional 4th column will be loaded as well and stored in memory too,
     * as the 'extra' attribute on the GRange objects.
     * <p>
     * Note: this method returns the same BED object as a result to enable chaining
     */
    @CompileStatic
    BED load(Map options=[:]) {
        
        if(options.containsKey('withExtra'))
            this.withExtra = options.withExtra
        
        if(!withExtra) {
          eachRange(options) { String chr, int start, int end ->
              add(chr,start,end)
          }
        }
        else {
          eachRange(options) { String chr, int start, int end, String extra ->
              add(chr,start,end,extra)
          }
        }
        isLoaded = true
        return this
    }
    
    /**
     * Add the specified range to this BED file.
     * <p>
     * NOTE: the 'end' is treated as exclusive of the range covered. This
     * is consistent with BED file notation, but different to Groovy ranges.
     * Since most methods in this class use Groovy ranges, you will generally
     * get a Range object out that has end one less than the end you put in
     * with add. Thus if you are iterating through one BED file and adding the
     * ranges to another with this method, you <i>must</i> add one to the
     * end position that you pass to this method!
     * <p>
     * Will trigger the flag to indicate this BED is loaded in memory.
     */
    @CompileStatic
    BED add(String chr, int start, int end, Object extra = null) {
		super.addRegion(chr,start,end,extra)
        isLoaded = true
        return this
    }
    
    /**
     * Iterate every single base in the regions of the specified bed file
     * @param fileName
     */
    static void eachPosition(String fileName, Closure c) {
        new BED(fileName).eachPosition(c)
    }
    
    static BED eachInputRange(Closure c) {
        def b = new BED(System.in)
        b.eachRange(c)
        return b
    }
    
    static BED parse() {
        def b = new BED(System.in)
        b.load()
        return b
    }
    
    static BED parse(String fileName) {
        BED b = new BED(fileName)
        b.load()
        return b
    }
	
	boolean isCase(IRegion region) {
		return this.getOverlaps(region.chr, region.range.from, region.range.to)
	}
   
    String toString() {
        int numRanges = getNumberOfRanges()
        if(numRanges>0) {
            "${numRanges} regions starting at ${this.iterator().next()}"
        }
        else {
            "Empty BED (no regions)"
        }
    }
}
