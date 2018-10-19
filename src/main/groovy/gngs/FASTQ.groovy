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

import java.nio.file.Files
import java.text.ParseException;

import java.util.zip.GZIPInputStream;

import groovy.transform.CompileStatic;
import groovy.transform.stc.ClosureParams
import groovy.transform.stc.SimpleType

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
    static void eachRead(@ClosureParams(value=SimpleType,options=['gngs.FASTQRead']) Closure c) {
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
    static void filter(String fileName1, String fileName2, Writer output, @ClosureParams(value=SimpleType,options=['gngs.FASTQRead']) Closure<Boolean> c) {
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
     * Filter R1 reads with from fileName1 and write them to 
     * the output in FASTQ format trimming the Chromium barcode off the front.
     *
     * @param fileName1
     * @param output
     * @param c
     */
    @CompileStatic
    static void filter10Xnobarcode(String fileName1, Writer output, Closure c) {
        int passed = 0
        int processed = 0
        FASTQ.eachRead(fileName1) { FASTQRead r1 ->
            ++processed
            def barcode10X = r1.trimEnd(r1.bases.size() - 16)
            def r1Rest = r1.trimEnd(0, 16)
            def result = c(r1Rest,barcode10X)
            if(result == true) {
                output.println(r1.header)
                output.println(r1Rest.bases)
                output.println("+")
                output.println(r1Rest.quals)
                ++passed
            }
        }
    }


    /**
     * Filter paired reads with index from fileName{1,2,3} and write them to 
     * the output in 10X (Chromium) Clariat compatible interleaved format.
     * 
     * @param fileName1
     * @param fileName2
     * @param fileName3
     * @param output
     * @param c
     */
    @CompileStatic
    static void filter10X(String fileName1, String fileName2, String fileName3, Writer output, Closure c) {
        int passed = 0
        int processed = 0
        FASTQ.eachPairWithIndex(fileName1,fileName2,fileName3) { FASTQRead r1, FASTQRead r2, FASTQRead i ->
            ++processed
            def barcode10X = r1.trimEnd(r1.bases.size() - 16)
            def r1Rest = r1.trimEnd(0, 16)
            def result = c(r1Rest,r2,barcode10X,i)
            if(result == true) {
                output.println(r1.header)
                output.println(r1Rest.bases)
                output.println(r1Rest.quals)
                output.println(r2.bases)
                output.println(r2.quals)
                output.println(barcode10X.bases)
                output.println(barcode10X.quals)
                output.println(i.bases)
                output.println(i.quals)
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
    static void filter(String fileName1, String fileName2, String output1, String output2, @ClosureParams(value=SimpleType,options=['gngs.FASTQRead','gngs.FASTQRead']) Closure c) {
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
    static void filter(String fileName1, String fileName2, Writer output1, Writer output2, @ClosureParams(value=SimpleType,options=['gngs.FASTQRead','gngs.FASTQRead']) Closure c) {
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
    static void eachPair(String fileName1, String fileName2, @ClosureParams(value=SimpleType,options=['gngs.FASTQRead','gngs.FASTQRead']) Closure c) {
        
        Utils.reader(fileName1) { Reader reader1 ->
          Utils.reader(fileName2) { Reader reader2 ->
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
    static void eachPairWithIndex(String fileName1, String fileName2, String fileName3, Closure c) {

        Utils.reader(fileName1) { Reader reader1 ->
          Utils.reader(fileName2) { Reader reader2 ->
            Utils.reader(fileName3) { Reader reader3 ->
              ProgressCounter counter = new ProgressCounter(withRate:true, withTime:true)
              try {
                  while(true) {
                    FASTQRead read1 = consumeRead(reader1)
                    FASTQRead read2 = consumeRead(reader2)
                    FASTQRead index = consumeRead(reader3)
                    if(read1 == null) // end of file
                        break
                    if(read2 == null)
                        throw new IllegalStateException("Trailing reads found in $fileName2 that are not present in $fileName1")
                    if(index == null)
                        throw new IllegalStateException("Trailing reads found in $fileName2 that are not present in $fileName1")

                    if((read1.name != read2.name) || (read1.name != index.name))
                        throw new IllegalStateException("Read $read1.name from $fileName1 is not matched by read at same line in $fileName2 ($read2.name) and in $fileName3 ($index.name). Reads need to be in same order in all files.")
                    c(read1,read2,index)
                    counter.count()
                  }
              }
              catch(Abort a) {
                  // expected
              }
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
    static void eachRead(String fileName, @ClosureParams(value=SimpleType,options=['gngs.FASTQRead']) Closure c) {
        Utils.reader(fileName) { Reader reader ->
            iterateReader(reader,c)
        }
    }
    
    
    @CompileStatic
    private static void iterateReader(Reader reader, Closure c) {
        // Read the file 4 lines at a time
        while(true) {
          final FASTQRead read = consumeRead(reader)
          if(read == null)
              return
                  
          c(read)
        }
    }
}
