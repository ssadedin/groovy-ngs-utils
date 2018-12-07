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
package gngs

import java.text.ParseException;

import java.util.zip.GZIPInputStream;

import groovy.transform.CompileStatic;

/**
 * Models a single FASTQ read, including the header (with name), bases and 
 * quality information.
 */
@CompileStatic
class FASTQRead {
    
    /**
     * This constructor is really just for testing purposes.
     * However you can use it if you want to fake all the data except the bases.
     * 
     * @param bases
     */
    FASTQRead(String bases) {
        this.name = Hash.sha1(bases)
        this.bases = bases
        this.quals = 'A' * bases.size()
        this.header = this.name
    }
    
    FASTQRead(String header, String bases, String quals) {
        this.header = header
        
        int slashIndex = header.indexOf('/')
        int spaceIndex = header.indexOf(' ')
        if(spaceIndex<0) {
            this.name = slashIndex>0?header.subSequence(0, slashIndex):header
        }
        else // There is a space, but if a slash comes before use that
        if(slashIndex>=0) {
            this.name = header.subSequence(0, Math.min(spaceIndex,slashIndex))
        }
        else
          this.name = header.subSequence(0, spaceIndex)
            
        this.bases = bases
        this.quals = quals
    }
    
    FASTQRead trimEnd(int count, int countStart=0) {
        new FASTQRead(header, bases.substring(countStart,bases.size()-count), quals.substring(countStart,quals.size() - count))
    }
    
    void write(Writer w) {
        w.println(header)
        w.println(bases)
        w.println("+")
        w.println(quals)
    }
    
    int size() {
        bases.size()
    }
    
    String header
    CharSequence name
    String bases
    String quals
}

/**
 * Simple support for reading fastq files.
 * <p>
 * This class supports reading either single end reads (via {@link #eachRead(Closure)}) or
 * paired end reads (via {@link #eachPair(String, String, Closure)}.
 * <p>
 * Gzipped files are supported transparently - if the file name ends in 'gz' then 
 * the file is assumed to be gzipped and decoded as such.
 * <p>
 * This class only supports reading FASTQ reads, not writing. However the individual reads
 * can be written via the {@link FASTQRead} class.
 * 
 * @author Simon
 */
class FASTQ {
    
    final static long PRINT_INTERVAL_MS = 15000
    
    @CompileStatic
    static void eachRead(Closure c) {
        // Cheat, fails on windows
        eachRead("/dev/stdin",c)
    }
    
    /**
     * Filter paired reads from fileName1 and fileName2 and write them to 
     * uncompressed output files with extensions .filter.fastq based on the 
     * input file names, only where true is returned from the given closure
     * 
     * @param fileName1
     * @param fileName2
     * @param c
     */
    @CompileStatic
    static void filter(String fileName1, String fileName2, Closure c) {
        new File(fileName1.replaceAll(/\.fastq.gz/, ".filter.fastq")).withWriter { Writer w1 ->
            new File(fileName2.replaceAll(/\.fastq.gz/, ".filter.fastq")).withWriter { Writer w2 ->
                filter(fileName1, fileName2, w1, w2, c)
            }
        }
        
    }
    
    /**
     * Filter paired reads from fileName1 and fileName2 and write them to 
     * the output in interleaved format.
     * 
     * @param fileName1
     * @param fileName2
     * @param c
     */
    @CompileStatic
    static void filter(String fileName1, String fileName2, Writer output, Closure c) {
        int passed = 0
        int processed = 0
        FASTQ.eachPair(fileName1,fileName2) { FASTQRead r1, FASTQRead r2 ->
            ++processed
            def result = c(r1,r2)
            if(result == true) {
                r1.write(output)
                r2.write(output)
                ++passed
            }
        }
    }
  
    /**
     * Filter paired reads from fileName1 and fileName2 and write them to 
     * uncompressed output files with extensions .filter.fastq based on the 
     * input file names, only where true is returned from the given closure
     * 
     * @param fileName1
     * @param fileName2
     * @param c
     */
    @CompileStatic
    static void filter(String fileName1, String fileName2, String output1, String output2, Closure c) {
        Utils.outputWriter(output1).withWriter { w1 ->
            Utils.outputWriter(output2).withWriter { w2 -> 
                filter(fileName1, fileName2, w1, w2, c)
            }
        }
    }
     
    /**
     * Filter paired reads from fileName1 and fileName2 and write them to 
     * uncompressed output files with extensions .filter.fastq based on the 
     * input file names, only where true is returned from the given closure
     * 
     * @param fileName1
     * @param fileName2
     * @param c
     */
    @CompileStatic
    static void filter(String fileName1, String fileName2, Writer output1, Writer output2, Closure c) {
        int passed = 0
        int processed = 0
        FASTQ.eachPair(fileName1,fileName2) { FASTQRead r1, FASTQRead r2 ->
            ++processed
            def result = c(r1,r2)
            if(result == true) {
                r1.write(output1)
                r2.write(output2)
                ++passed
            }
        }
    }
  
    
    @CompileStatic
    static void eachPair(String fileName1, String fileName2, Closure c) {
        
        openStream(fileName1).withReader { Reader reader1 ->
          openStream(fileName2).withReader { Reader reader2 ->
            ProgressCounter counter = new ProgressCounter(withRate:true, withTime:true)
            try {
                while(true) {
                  FASTQRead read1 = consumeRead(reader1)
                  FASTQRead read2 = consumeRead(reader2)
                  if(read1 == null) // end of file
                    break
                  if(read2 == null) 
                      throw new IllegalStateException("Trailing reads found in $fileName2 that are not present in $fileName1")
                      
                  if(read1.name != read2.name)
                      throw new IllegalStateException("Read $read1.name from $fileName1 is not matched by read at same line in $fileName2 ($read2.name). Reads need to be in same order in both files.")
                  c(read1,read2)
                  counter.count()
                }
            }
            catch(Abort a) {
                // expected
            }
          }
        }
    }
    
    @CompileStatic
    static FASTQRead consumeRead(Reader reader) {
          String name = reader.readLine()
          if(name == null)
              return
          String bases = reader.readLine()
          if(bases == null)
              throw new ParseException("Incorrect FASTQ format: no bases after read name $name", -1)
              
          String sep = reader.readLine()
          if(sep != "+")
              throw new ParseException("Incorrect FASTQ format: expected '+' on line after read name $name", -1)
          String quals = reader.readLine()
          if(quals == null)
              throw new ParseException("Incorrect FASTQ format: no quality socres after read name $name", -1)
              
          return new FASTQRead(name, bases, quals)
    }
    
    @CompileStatic
    static void eachRead(String fileName, Closure c) {
        
        InputStream readIn = openStream(fileName)
            
        readIn.withReader { Reader reader ->
            // Read the file 4 lines at a time
            while(true) {
              FASTQRead read = consumeRead(reader)
              if(read == null)
                  break
                  
              c(read)
            }
        }
    }
    
    static InputStream openStream(String fileName) {
        int bufferSize = 1024*1024
        if(fileName.endsWith(".gz"))
          new BufferedInputStream(new GZIPInputStream(new FileInputStream(fileName)), bufferSize)
        else
          new BufferedInputStream(new File(fileName).newInputStream(), bufferSize)
    }
}
