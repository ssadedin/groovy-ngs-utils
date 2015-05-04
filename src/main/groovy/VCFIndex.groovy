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
 * Support for querying indexed VCF files via random access.
 * <p>
 * In general, if you only want to access a small region of a VCF (eg: to look up variants in 
 * a specific region, etc) then this class is your best bet. If, alternatively, you have a 
 * smaller VCF file (say, less than 100,000 lines) and you expect to query a significant number 
 * of the variants in it, you should go straight for the {@link VCF} class instead and just
 * load everything into memory.
 * <p>
 * The VCFIndex class is highly optimised for random access to VCF contents. It 
 * uses memory mapped files to buffer the VCF to memory to minimize the actual hits to the 
 * file system.
 * <p>Using the VCFIndex class to look up variants in a region is very simple:
 * <pre>
 * VCFIndex index = new VCFIndex("test.vcf")
 * index.query("chrX",1000,2000) { Variant v ->
 *  println "Variant $v is in the range chrX:1000-2000"
 * }
 * </pre>
 * You can also easily test if a Variant exists within the VCF file:
 * <pre>
 * Variant v = new Variant(chr:"chrX", alt:"A", ref:"T")
 * VCFIndex index = new VCFIndex("test.vcf")
 * if(v in index)
 *      println "Variant $v is found in test.vcf"
 * </pre>
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class VCFIndex {
    
    /**
     * File name of VCF file
     */
    String fileName
    
    /**
     * Raw file corresponding to the VCF
     */
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
    
    /**
     * Find the first variant in an interval that returns true from the
     * given closure
     * <p>Example:
     * <pre>
     * VCFIndex index = new VCFIndex("test.vcf")
     * Variant v = index.find("chrX", 1000,2000) { it.type == "DEL" && it.size() > 4 }
     * println "First deletion greater than 4 bases in size is $v"
     * </pre>
     * 
     * @param chr   Chromosome to scan
     * @param start Starting position in chromosome
     * @param end   End position in chromosome
     * @param c     Closure to call, passing each Variant from the VCF as the first argument.
     *              
     * @return      The first Variant in the specified range for which the Closure c returns an
     *              argument that does not evaluate to false (null, 0, etc.)
     */
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
    
    boolean isCase(Variant v) {
        this.contains(v)
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
    
    void query(IRegion r, Closure c) {
        query(r.chr, r.from, r.to, c)
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
            this.vcfBuffers = (List<MappedByteBuffer>)[0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L].collect { long i ->
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
    
    
    /**
     * Attempts to locate the given Annovar variant in this VCF file.
     * <p>
     * Annovar is a popular annotation tool which outputs results in its own 
     * proprietary CSV format. Unfortunately, this format does not use the same 
     * reference coordinates, and thus an Annovar variant is not directly translatable 
     * to a VCF variant, and it is necessary to scan a range of indexes in an original VCF file
     * to locate the corresponding Annovar variant. This function implements the required logic.
     * 
     * @param chr           chromosome / reference sequence
     * @param start         Annovar variant start position
     * @param end           Annovar variant end position
     * @param windowSize    the size of window to scan. VCF files can include overlapping variants in 
     *                      arbitrarily long windows. Thus it is necessary to start scanning the VCF
     *                      significantly before the location of the actual start of the DNA change
     *                      to be sure of finding the Annovar variant. In general, the window size
     *                      should be as large as the largest indels you expect to have in your
     *                      VCF file.
     * 
     * @return  if the variant is found, a map containing "variant" and "allele" entries.
     *          the variant is the Variant object, while the allele is the zero-based index of 
     *          the allele within the  Variant that corresponds to the Annovar variant specified.
     */
    Map findAnnovarVariant(String chr, def start, def end, String obs, int windowSize=15) {
        start = start.toInteger()
        end = end.toInteger()
        int alleleIndex = -1;
        Variant v = this.find(chr,start-windowSize,end) { variant -> alleleIndex = variant.equalsAnnovar(chr,start,obs) }
        return v ? [variant:v, allele: alleleIndex-1] : null
    }
}
