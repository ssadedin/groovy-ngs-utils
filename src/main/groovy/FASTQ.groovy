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

import java.text.ParseException;
import java.util.zip.GZIPInputStream;

import groovy.transform.CompileStatic;

@CompileStatic
class FASTQRead {
    FASTQRead(String header, String bases, String quals) {
        this.header = header
        
        int slashIndex = header.indexOf('/')
        int spaceIndex = header.indexOf(' ')
        if(spaceIndex<0) {
            this.name = slashIndex>0?header.subSequence(0, slashIndex):header
        }
        else
            this.name = header.subSequence(0, Math.min(spaceIndex,slashIndex))
            
        this.bases = bases
        this.quals = quals
    }
    
    FASTQRead trimEnd(int count) {
        new FASTQRead(header, bases.substring(0,bases.size()-count), quals.substring(0,quals.size() - count))
    }
    
    void write(Writer w) {
        w.println(header)
        w.println(bases)
        w.println("+")
        w.println(quals)
    }
    
    String header
    CharSequence name
    String bases
    String quals
}

class FASTQ {
    
    final static long PRINT_INTERVAL_MS = 15000
    
    @CompileStatic
    static void eachRead(Closure c) {
        // Cheat, fails on windows
        eachRead("/dev/stdin",c)
    }
    
    @CompileStatic
    static void eachPair(String fileName1, String fileName2, Closure c) {
        
        openStream(fileName1).withReader { Reader reader1 ->
          openStream(fileName2).withReader { Reader reader2 ->
            def progress =new ProgressCounter()
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
              progress.count()
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
            
        def progress =new ProgressCounter()
        readIn.withReader { Reader reader ->
            // Read the file 4 lines at a time
            while(true) {
              FASTQRead read = consumeRead(reader)
//              c(read.name,read.bases,read.quals)
              if(read == null)
                  break
                  
              c(read)
              progress.count()
            }
        }
    }
    
    static InputStream openStream(String fileName) {
        if(fileName.endsWith(".gz"))
          new GZIPInputStream(new FileInputStream(fileName))
        else
          new File(fileName).newInputStream()
    }
}
