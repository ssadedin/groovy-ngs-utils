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

import groovy.transform.CompileStatic;


/**
 * Simple utility class to assist in parsing BED files
 * <p>
 * At the moment it is limited (mainly) to working with human data because 
 * chromosomes chrX, chrY and chrM are handled as special cases. Also, 
 * other contigs and random chromosomes are not handled either.
 * <p>
 * Some methods hit the file directly, while others rely on reading it 
 * into memory. If you don't care too much about performance you can 
 * not care aboaut this. But if you care about performance you should 
 * check how the method is going to operate and choose accordingly.
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class BED {
    
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
     * Index for looking up overlaps
     */
    Map<String, RangeIndex> index = [:]
    
    /**
     * A list of ranges in the order they were loaded
     */
    Map<String,List> allRanges = [:]
    
    /**
     * The stream from which bed line will be read
     */
    InputStream bedFileStream = null
    
    /**
     * File from which bed is read
     */
    File bedFile
    
    /**
     * Map containing all the ranges loaded from the bed file,
     * index by start position, and partially reduced to make
     * lookups efficient.
     * Note: only loaded if load() is called!
     */
    TreeMap<Long,Range> startRanges = new TreeMap()
    
    /**
     * Map containing all the ranges loaded from the bed file,
     * indexed by end position, and partially reduced to make
     * lookups efficient.
     */
    TreeMap<Long,Range> endRanges = new TreeMap()
    
    /**
     * Empty bed file
     */
    BED(Map attributes=[:]) {
    }
    
    BED(Map attributes=[:], String fileName) {
        this(new File(fileName))
    }
    
    BED(Map attributes=[:], File file, Closure c = null) {
        this.bedFile = file
        this.bedFileStream = new FileInputStream(file)
        if(attributes && attributes.withExtra != null)
            this.withExtra = attributes.withExtra
    }
    
    BED(InputStream inStream, Closure c = null) {
        this.bedFileStream = inStream
    }
    
    @CompileStatic
    int size() {
        int sizeCount = 0
        this.eachRange { String chr, int start, int end ->
            sizeCount += (end - start + 1)
        }
        return sizeCount
    }
    
    /**
     * Iterate over each line in the BED file and invoke the specified
     * closure with bed file information for the line
     */
    @CompileStatic
    void eachRange(Closure c) {
        
        if(bedFileStream == null && isLoaded) {
            eachLoadedRange(c)
            return
        }
        
        int count = 0
        def includeInfo = c.getMaximumNumberOfParameters()>3;
        
        try {
            if(bedFileStream == null || bedFileStream.available()<0) {
                if(bedFile != null) {
                  bedFileStream = new FileInputStream(bedFile)
                }
                
            }
        }
        catch(IOException ex) {
            try { bedFileStream.reset() } catch(Exception e) {}
            if(bedFile != null) 
                  bedFileStream = new FileInputStream(bedFile)
         }
        
        if(bedFileStream == null)
            return
        
        def abort = { throw new Abort("break from BED file iteration") }
        ProgressCounter progress = new ProgressCounter(lineInterval:10, timeInterval:10000)
        try {
          bedFileStream.eachLine { String line ->
              ++count
              if(line.startsWith('#'))
                  return
              if(line.startsWith("track"))    
                  return
              if(line.startsWith("browser"))    
                  return
              String [] fields = line.split('\t')
              String chr = fields[0]
              int start = -1, end = -1
              try {
                start = Integer.parseInt(fields[1])
                end = Integer.parseInt(fields[2])
              }
              catch(Exception e) {
                  System.err.println "Unable to parse line $count : \n$line\n\nException reported: $e"   
                  throw e
              }
              try {
                if(includeInfo==true)
                    c.call(chr,start,end,fields[3])
                else
                    c.call(chr,start,end)
                    
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
    }
    
//    @CompileStatic
    void eachLoadedRange(Closure c) {
        boolean includeInfo = c.getMaximumNumberOfParameters()>3;
        this.allRanges.each { String chr, List ranges ->
            for(groovy.lang.IntRange r in ranges) {
              if(includeInfo) {
                  c.call(chr,r.from,r.to,r instanceof GRange ? r.extra : null)
              }
              else
                  c.call(chr,r.from,r.to)
            }
        }
    }
    
    List<Range> intersect(String chr, int start, int end) {
        RangeIndex chrIndex = this.index[chr]
        if(chrIndex == null)
            return []
        
        return chrIndex.intersect(start,end)
    }
    
    List<Range> getOverlaps(String chr, int start, int end) {
        RangeIndex chrIndex = this.index[chr]
        if(chrIndex == null)
            return []
        return chrIndex.getOverlaps(start,end)
    }
    
    List<Range> subtractFrom(String chr, int start, int end) {
        RangeIndex chrIndex = this.index[chr]
        if(chrIndex == null)
            return [start..end-1]
        
        return chrIndex.subtractFrom(start,end)
    }
    
    /**
     * Call closure c for each range that overlaps the given position
     */
    @CompileStatic
    void eachOverlap(String chr, int pos, Closure c) {
        if(!isLoaded) {
          eachRange { String rangeChr, int rangeStart, int rangeEnd ->
              if(rangeChr == chr && pos >= rangeStart && pos << rangeEnd) {
                  c(rangeChr, rangeStart, rangeEnd)
              }
          }
        }
        else {
            RangeIndex chrIndex = this.index[chr]
            if(chrIndex == null)
                return
            for(groovy.lang.Range r in chrIndex.getOverlaps(pos)) {
                c(chr, r.from, r.to)
            }
        }
    }
    
    List<Range> startingAt(String chr, int pos) {
        RangeIndex chrIndex = this.index[chr]
        if(chrIndex == null)
            return []
        
        return chrIndex.startingAt(pos)
    }
    
    List<Range> endingAt(String chr, int pos) {
        RangeIndex chrIndex = this.index[chr]
        if(chrIndex == null)
            return []
        
        return chrIndex.endingAt(pos)
    }
    
    /**
     * Returns the prior range that has its end closest
     * to the given position.
     */
    groovy.lang.IntRange previousRange(String chr,int pos) {
        RangeIndex chrIndex = this.index[chr]
        if(chrIndex == null)
            return null
        return chrIndex.previousRange(pos)
    } 
    
    /**
     * Returns the next range that has its beginning
     * closest to the given position.
     */
    groovy.lang.IntRange nextRange(String chr, int pos) {
        RangeIndex chrIndex = this.index[chr]
        if(chrIndex == null)
            return null
        return chrIndex.nextRange(pos)
    }
    
    @CompileStatic
    void remove(String chr, Range r) {
        RangeIndex chrIndex = this.index[chr]
        if(!chrIndex) 
            throw new IllegalArgumentException("Cannot remove region from chrosomome $chr : it does not have any regions")
        chrIndex.remove(r)
        
        if(this.allRanges[chr]) // Could be really slow?
            this.allRanges[chr].removeAll { Range r2 -> r.from == r2.from && r.to == r2.to }
    }
    
    /**
     * Returns the range that is "closest" to the given
     * position.
     * 
     *   1. if no range overlaps, the closer of a) the end
     *      of the nearest prior range vs b) the start of the
     *      nearest following range is returned. 
     *      
     *   2. if one or more ranges overlap, returns the overlapping
     *      range with the start or end closest to the specified position
     *   
     * @param chr
     * @param pos
     * @return
     */
    groovy.lang.IntRange nearest(String chr, int pos) {
        this.index[chr]?.nearest(pos)
    }
    
    List<Object> getExtrasAtPosition(String chr, int position) {
        List result = []
        
        RangeIndex chrIndex = this.index[chr]
        if(chrIndex == null)
            return result
        for(groovy.lang.Range range in chrIndex.getOverlaps(position)) {
            if(range instanceof GRange)
                result << range.extra
        } 
        return result
    }
    
    /**
     * Iterate each position in the BED file and call the given closure with it
     * 
     * @TODO: merge overlapping regions
     */
    @CompileStatic 
    void eachPosition(Closure c) {
        int totalSize = size()
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
    BED load() {
        if(!withExtra) {
          eachRange { String chr, int start, int end ->
              add(chr,start,end)
          }
        }
        else {
          eachRange { String chr, int start, int end, String extra ->
              add(chr,start,end,extra)
          }
        }
        isLoaded = true
        return this
    }
    
    /**
     * Add the specified range to this BED file.
     * <p>
     * NOTE: the 'end' is treated as exclusive of the range covered.
     * <p>
     * Will trigger the flag to indicate this BED is loaded in memory.
     */
    @CompileStatic
    BED add(String chr, int start, int end, Object extra = null) {
        
        // Add to full range index
        RangeIndex chrIndex = index[chr]
        if(!chrIndex) {
            chrIndex = new RangeIndex()
            index[chr] = chrIndex
        }
        
        groovy.lang.IntRange newRange = null
        if(extra == null && !withExtra) {
           newRange = start..(end-1)
        }
        else {
           newRange = new GRange(start,end-1,extra)   
        }
        chrIndex.add(newRange)
        
        if(!allRanges.containsKey(chr))
            allRanges[chr] = new ArrayList(1000)
            
        allRanges[chr] << newRange
        
        // Add to reduced index for fast lookups 
        // TODO: remove
        long base = chrToBaseOffset(chr)
        GRange r = new GRange(start, end, extra)

        long startKey = base + start
        if(startRanges.containsKey(startKey)) {
            if(r.to < end) {
                assert false : "Should never come here!"
                r = new GRange(r.from,end,r.extra)
            }
        }
        startRanges[startKey] = r
        endRanges[base + end] = r
        isLoaded = true
        return this
    }
    
    @CompileStatic
    private long chrToBaseOffset(String chr) {
        if(chr == "chrX") {
            return (23L << 32)
        }
        else
        if(chr == "chrY") {
            return (24L << 32)
        }
        else
        if(chr == "chrM") {
            return (25L << 32)
        }
        else
        if(chr.startsWith("chr")) {
            return ((chr.substring(3).toLong()-1) << 32)
        }
        else {
            return ((chr.toLong()-1) << 32)
        }
    }
    
    /**
     * Return true iff the given chromosome and position fall into 
     * a range covered by this BED file
     */
    boolean contains(String chr, int position) {
        RangeIndex chrIndex = index[chr]
        return chrIndex != null && (position in chrIndex)
    }
    
    boolean isEmpty() {
        return allRanges.isEmpty()
    }
    
    /**
     * Iterate every single base in the regions of the specified bed file
     * @param fileName
     */
    static void eachPosition(String fileName, Closure c) {
        new BED(fileName).eachPosition(c)
    }
    
    static void each(Closure c) {
        new BED(System.in).eachRange(c)
    }
    
    static BED parse(String fileName) {
        BED b = new BED(fileName)
        b.load()
        return b
    }
}
