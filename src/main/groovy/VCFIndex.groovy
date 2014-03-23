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

import java.nio.MappedByteBuffer
import java.nio.channels.FileChannel.MapMode

import org.broad.tribble.index.Block
import org.broad.tribble.index.Index
import org.broad.tribble.index.IndexFactory
import org.codehaus.groovy.runtime.StackTraceUtils


/**
 * Helper for reading and querying VCF files that are indexed
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class VCFIndex {
    
    String fileName
    
    RandomAccessFile indexedFile
    
    /**
     * A 256 kb memory buffer that will try to minimize operations
     * on the indexed file
     */
    public static final int BUFFER_SIZE = 1024 * 256;
    
    public static final int ONE_GIG = 1024*1024*1024
    
    List<MappedByteBuffer> vcfBuffers = null
    
    /**
     * Indexes used for random access to large VCFs. Used with *.query() functions
     */
    Index index = null
    
    /**
     * A dummy VCF used to parse / hold the header information
     */
    VCF headerVCF = null

    /**
     * Load and index for the given VCF file
     * 
     * @param fileName
     */
    public VCFIndex(String fileName) {
         this.fileName = fileName
         
         if(!new File(fileName).exists()) 
             throw new FileNotFoundException(new File(fileName).absolutePath, "The VCF file specified does not exist or could not be accessed")
         
         File indexFile = new File(fileName + ".idx")
         if(!indexFile.exists()) 
             throw new FileNotFoundException(indexFile.absolutePath, 
                 """
                     The VCF file specified is not indexed. Please index the VCF file, for example, using: 
                 
                     igvtools index $fileName
                 
                 """.stripIndent())
             
         this.index = IndexFactory.loadIndex(fileName + ".idx")
         this.headerVCF = new VCF(fileName)
     }
    
    Variant find(String chr, int start, int end, Closure c) {
        Variant result = null
        query(chr,start,end) { Variant v ->
            if(c(v)) {
                result = v
                return false
            }
        }
        return result
    }
    
    /**
     * Returns true if this VCF contains the specified variant
     * <p>
     * IMPORTANT: the other variant is required to contain only
     *            a single allele for this function to return correct 
     *            results!
     * 
     * @param v Variant to test for
     * @return  The variant iff this VCF contains the variant, null otherwise
     */
    Variant contains(Variant v) {
        Variant found = this.find(v.chr, v.pos-10, v.pos + v.size()+10) { myv ->
            myv.pos == v.pos && myv.allelesAndTypes.any { it[0] == v.alt && it[1] == v.type }
        }
        
        return found
    }
    
    /**
     * Query the current VCF file for all variants in the specified
     * region
     *
     * @param chr   Chromosome / Sequence name
     * @param start Start position
     * @param end   End position
     * @param c     Closure to call back
     */
    @CompileStatic
    void query(String chr, int start, int end, Closure c) {
        
        if(this.index == null) {
            throw new RuntimeException("Random access queries require a VCF file loaded into memory")
        }
        else {
          if(this.indexedFile == null) {
            this.indexedFile = new RandomAccessFile(new File(this.fileName), "r")
            
            // Map the whole damn thing into memory, with overlapping buffers
            this.vcfBuffers = (List<MappedByteBuffer>)[0L, 1L, 2L, 3L, 4L].collect { long i ->
                i*ONE_GIG < indexedFile.length() ? indexedFile.channel.map(MapMode.READ_ONLY, i*ONE_GIG, Math.min(indexedFile.length() - i*ONE_GIG,(int)(1.5*ONE_GIG))) : null
            }
          }
            
          List blocks = null
          try {
            blocks = this.index.getBlocks(chr, start, end)
          }
          catch(IllegalArgumentException ex) {
             // This occurs when there are none of the requested chromosome in the VCF file -
             // Just treat it as if nothing was found
              return
          }
          
          for(Block block in blocks) {
              
              byte [] buffer = new byte[block.size+1]
              
              int bufferIndex = (int)Math.floor(block.startPosition / (long)ONE_GIG)
              
              if(bufferIndex > vcfBuffers.size())
                  throw new RuntimeException("Trying to access position beyond end of buffer $block.startPosition maps to buffer $bufferIndex")
              
              if(bufferIndex<0)
                  throw new RuntimeException("Trying to access negative position?")
              
//              System.err.println "Trying for start position $block.startPosition, index = $bufferIndex"
              
              MappedByteBuffer vcfBuffer =  (MappedByteBuffer)vcfBuffers.get(bufferIndex)
              vcfBuffer.position((int)(block.startPosition - bufferIndex*ONE_GIG))
              vcfBuffer.get(buffer, 0, (int)(block.size))
              
              boolean cont = new ByteArrayInputStream(buffer).withReader("ASCII") { Reader r ->
                String line
                int blockLineCount = 0
                while((line = r.readLine()) != null) {
//                    println "Parsing line: $line"
                    
                    line = line.trim()
                    
                    if(line.startsWith('#'))
                        continue
                        
                    if(line.isEmpty())
                        break
                        
                    ++blockLineCount
                    try {
                        Variant v = Variant.parse(line)
                        if(v.pos > end)
                            break

                        if(v.pos < start)
                            continue

                        v.header = this.headerVCF
                        try {
                          if(c(v)==false)
                              return false
                        }
                        catch(Exception e) {
                          System.err.println("Failure in processing line $blockLineCount [$line] empty=${line.trim().isEmpty()} lineSize=${line.size()} start position $block.startPosition size=$block.size buffer index $bufferIndex")
//                          e.printStackTrace()
                          StackTraceUtils.deepSanitize(e).printStackTrace()
//                          StackTraceUtils.printSanitizedStackTrace(e)
                        }
                    }
                    catch(Exception e) {
                        System.err.println("Failure in parsing line $blockLineCount [$line] empty=${line.trim().isEmpty()} lineSize=${line.size()} start position $block.startPosition size=$block.size buffer index $bufferIndex")
                    }
                }
                return true;
            }
            if(!cont)
                break
          }
        }
    }
}
