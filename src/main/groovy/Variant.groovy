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

import java.nio.ByteBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel.MapMode;

import groovy.transform.CompileStatic;

import org.broad.tribble.index.Block;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.codehaus.groovy.runtime.StackTraceUtils;

/**
 * Support for parsing annotations from SnpEff
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class SnpEffInfo {
    String type
    String impact
    String gene
    String transcript
    
    // The original subclause within the info line from which
    // this SnpEffInfo object was parsed
    String info
    
    String toString() {
        "type=$type,gene=$gene,rank=$impact"
    }
}

/**
 * The VCF class supports simple parsing, filtering, querying and updating of 
 * VCF files.
 * <p>
 * The primary method is the {@link #parse(Closure)} method which can be 
 * optionally passed a file name, stream, or if left out, will parse
 * from standard input. When a Closure is supplied, it can be used for
 * filtering and updating the VCF file as it is read. Each record is passed
 * to the closure as it is read, and if false is returned
 * from the closure it will be discarded. Various forms of information can
 * be queried - the INFO section is parsed to a map and genotype information
 * is available as well as computed dosage for each variant (0,1,2). Support
 * for SnpEff annotations is also provided with a dedicated snpEffInfo property
 * that returns a list of {@link SnpEffInfo} objects describing the various
 * annotations available.
 * <p>
 * Most information is lazily computed so the cost of parsing complex fields is 
 * deferred unless they are asked for.
 * <p>
 * VCF rows can be updated by entering an update closure:
 * <code>v.update { v.info.FOO='bar' }</code>
 * <p>
 * Limited support for "Pedigrees" is available. At the moment only grouping of
 * samples into families is supported and intra-family relationships are not modeled.
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class VCF {
    
    List<String> headerLines = []
    
    String [] lastHeaderLine
    
    List<String> samples 
    
    List<Variant> variants
    
    List<Pedigree> pedigrees
    
    List<Pedigree> samplePedigrees 
    
    String fileName
    
    RandomAccessFile indexedFile
    
    /**
     * Indexes used for random access to large VCFs. Used with *.query() functions
     */
    Index index = null
    
    VCF() {
    }
    
    VCF(String fileName) {
         this.index = IndexFactory.loadIndex(fileName + ".idx") 
         this.fileName = fileName
         
         new File(fileName).withReader { Reader r ->
           String line = r.readLine();
           while(line != null) {
             if(line.startsWith('#')) {
                 this.headerLines.add(line)
             }
             else
                 break
                 
             line = r.readLine()
           } 
         }
         parseLastHeaderLine()
    }
    
    /**
     * A 256 kb memory buffer that will try to minimize operations
     * on the indexed file
     */
    public static final int BUFFER_SIZE = 1024 * 256;
    
    public static final int ONE_GIG = 1024*1024*1024
    
    List<MappedByteBuffer> vcfBuffers = null
    
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

                        v.header = this
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
     * Convenience method to accept string for parsing file
     */
    static VCF parse(String fileName, List<Pedigree> peds = null, Closure c = null) {
        if(fileName == "-")
          parse(System.in,peds,c)
        else
          parse(new File(fileName),peds,c)
    }
    
    static VCF parse(String fileName, Closure c) {
        parse(fileName,null,c)
    }
    
    static VCF parse(File f, List<Pedigree> peds = null, Closure c = null) {
        parse(new BufferedInputStream(new FileInputStream(f)),peds,c)
    }
    
    static VCF parse(Closure c = null) {
        parse("-",null,c)
    }
    
    static VCF parse(List<Pedigree> peds, Closure c = null) {
        parse("-",peds,c)
    }
    
    static VCF parse(BufferedInputStream f, List<Pedigree> peds = null, Closure c = null) {
        parseImpl(f,false, peds,c)
    }
    
    static VCF filter(File f, Closure c = null) {
        filter(f,null,c)
    }
    
    static VCF filter(File f, List<Pedigree> peds, Closure c = null) {
        filter(new BufferedInputStream(new FileInputStream(f)),peds,c)
    }
    
    static VCF filter(String fileName, Closure c = null) {
        filter(fileName,null,c)   
    }
    
    static VCF filter(String fileName, List<Pedigree> peds, Closure c = null) {
        if(fileName == "-")
          filter(System.in,peds,c)
        else
          filter(new File(fileName),peds,c)
    }
    
    static VCF filter(Closure c = null) {
        filter("-",null,c)
    }
    
    static VCF filter(BufferedInputStream f, List<Pedigree> peds = null, Closure c = null) {
        parseImpl(f,true,peds,c)
    }
    
    /**
     * Parse the given VCF file, returning a full VCF, optionally filtered
     * by the closure (only rows returning true from closure will be included)
     * 
     * @param f file to parse
     * @param peds List of pedigree objects describing sample - family relationships
     * @param c filter closure
     */ 
    @CompileStatic
    private static VCF parseImpl(BufferedInputStream f, boolean filterMode, List<Pedigree> peds, Closure c) {
        VCF vcf = new VCF(pedigrees:peds)
        List<Variant> variants = []
        String lastHeaderLine
        long lastPrintTimeMs = System.currentTimeMillis()
        int count = 0
        boolean flushedHeader = false
        f.eachLine { String line ->
            
            ++count
            if(count % 10000 == 0 || System.currentTimeMillis() - lastPrintTimeMs > 30000) {
              System.err.println("Processed $count lines")
              lastPrintTimeMs = System.currentTimeMillis()
            }
            
            if(line.startsWith('#')) {
                vcf.headerLines.add(line)
                return
            }
            
            if(vcf.lastHeaderLine == null) {
                vcf.parseLastHeaderLine()
            }
            Variant v = Variant.parse(line)
            v.header = vcf
            try {
              if(!c || !(c(v)==false)) {
                  if(filterMode) {
                      if(!flushedHeader)  {
                          vcf.printHeader()
                          flushedHeader = true
                      }
                      println v.line
                  }
                  else
                    variants << v 
            }
          }
          catch(Exception e) {
             throw new RuntimeException("Procesing failed at line $count :\n" + line, e) 
          }
        }
        
        if(filterMode && !flushedHeader)  {
            vcf.printHeader()
            flushedHeader = true
        }
        
        vcf.variants = variants
        return vcf
    }
    
    Pedigree findPedigreeBySampleIndex(int i) {
        def peds = getPedigrees()
         if(samplePedigrees == null) {
             samplePedigrees = []
             samples.eachWithIndex { String s, int index ->
                 samplePedigrees[index] =  peds.find { p -> p.samples.contains(s) }
             }
         }    
         return samplePedigrees[i]
    }

    @CompileStatic
    int sampleIndex(String sampleName) {
        return samples.indexOf(sampleName)
    }
    
    /**
     * Extract sample names and other info from the final header line in the VCF file
     */
    void parseLastHeaderLine() {
        this.lastHeaderLine = this.headerLines[-1].split('[\t ]{1,}')
        this.samples = []
        for(int i=9; i<this.lastHeaderLine.size(); ++i) {
            this.samples.add(this.lastHeaderLine[i])
        }
    }
    
    /**
     * Index of variants by affected gene
     * Note that a variant can affect more than one gene, so
     * will appear multiple times in the value side of the map
     */
    Map<String,List<Variant>> genes = null
    
    @CompileStatic
    Map<String,List<Variant>> getGenes() {
        if(genes == null) {
          genes = [:]
          for(Variant v in variants) {
              for(SnpEffInfo eff in v.snpEffInfo.grep { SnpEffInfo inf -> inf.gene }) {
                List geneVars = genes[eff.gene]
                if(geneVars == null) {
                    genes[eff.gene] = []
                }
                genes[eff.gene] << v
              }
          }
        }
        return genes
    }
    
    void printHeader() {
        System.out.println(headerLines.join('\n'))
    }
    
    void print() {
        printHeader()
        for(Variant v in variants) {
            println(v.line)
        }
    }
}

/**
 * Help with parsing VCFs
 * @author simon.sadedin@mcri.edu.au
 */
class Variant {
    
    String chr
    int pos
    String ref
    String alt
    String [] alts
    String id
    String type
    String info
    float qual
    byte altByte
    VCF header

    Map<String,String> infos
    
    List<SnpEffInfo> snpEffInfo
    
    List<Integer> dosages
    
    Set<Pedigree> pedigrees
    
    List<Map<String,Object>> genoTypes
    
    String line
    
    String [] genoTypeFields
    
    // @CompileStatic
    void parseFields(String line) {
        
        def fields = line.split('[\t ]{1,}')
        
        id = fields[2]
        ref = fields[3]
        alts = fields[4].split(",")
        alt = alts[0]
        qual = (fields[5] == ".") ? -1.0f : Float.parseFloat(fields[5])
        chr = fields[0]
        info = fields[7]
        
        if(fields.length > 8) {
            
          if(genoTypeFields == null)
              genoTypeFields = fields[8].split(':')
          
          // Split the genoTypes field into separate values and parse them out
          genoTypes = fields[9..-1].collect { String gt -> 
              [genoTypeFields, gt.split(':')].transpose().collectEntries { it }
          } 
        }
          
        pos = Integer.parseInt(fields[1])
        
        type = "SNP"
        if(ref.size() < alt.size())
            type = 'INS'
        else
        if(ref.size() > alt.size())
            type = 'DEL'
        
         altByte = (byte)alt.charAt(0)
    }
    
    boolean snpEffDirty = false
    
    /**
     * Allows various fields to be updated and then synchronises the 
     * rest of the data with those updated fields
     * Not all fields are supported! see the ones that are set below.
     * <p>
     * The only update to snpEFF information is to remove individual 
     * annotations.
     */
    void update(Closure c) {
        this.snpEffDirty = false
        c()
        def fields = line.split('[\t ]{1,}')
        fields[0] = chr
        fields[1] = String.valueOf(pos)
        fields[3] = ref
        fields[4] = alt
        fields[5] = String.format("%2.2f",qual)
        
        if(snpEffDirty) {
            // The only operation we support on the snpEff meta data is to remove some.
            // Hence we just need to rebuild it by concatenating it
            getInfo().EFF = this.snpEffInfo*.info.join(",")
        }
        
        fields[7] = getInfo().collect { k,v -> "$k=$v"}.join(';')
        
        line = fields.join('\t')
        
        // Check if we need to add an INFO header line
        getInfo().collect { k,v ->
            
            if(!header.headerLines.find { it.startsWith("##INFO=<ID=$k") }) {
                
                int lastInfo = header.headerLines.findLastIndexOf { it.startsWith("##INFO=") }
                if(lastInfo < 0)
                    lastInfo = 1
                
                String valueType = "String"
                if(v instanceof Integer) {
                    valueType = "Integer"
                }
                else
                if(v instanceof Float) {
                    valueType = "Float"
                }
                
                header.headerLines = header.headerLines[0..lastInfo] + 
                    ["##INFO=<ID=$k,Number=1,Type=${valueType},Description=\"Added by groovy Variant processor\">"] +
                    header.headerLines[(lastInfo+1)..-1]
            }
        }
        this.@info = null
    }
    
    /**
     * VCF requires genotype to be the 1st field, but some tools (R, grrrr) write it in whatever order they 
     * feel like. This function corrects it.
     */
    void fixGenoTypeOrder() {
        def fields = line.split('[\t ]{1,}')
        
        def gtFields = fields[8]
        def gtInfo = fields[9]
        if(!gtFields.startsWith('GT:')) {
            def gtInfos = [gtFields.split(':'),gtInfo.split(':')].transpose().collectEntries { it }
            
            fields[8] = 'GT:' + gtInfos.collect { fieldName, fieldValue -> fieldName }.grep { it != 'GT' }.join(':')
            fields[9] = gtInfos['GT']+ ':' + gtInfos.grep { it.key != 'GT' }.collect { it.value }.join(':')
            line = fields.join('\t')
        }
    }
    
    @CompileStatic
    Map<String,String> getInfo() {
        if(infos == null) {
            infos =  [:]
            for(String s in info.split(';')) {
                int i = s.indexOf('=')
                if(i<0) {
                    infos[s]=Boolean.TRUE
                    continue
                }
                infos[s.substring(0,i)] = s.substring(i+1)
            }
        }
        return infos
    }
    
    /**
     * Return a list of SnpEffInfo objects, each describing
     * a separate SnpEff effect caused by the variant
     * @return
     */
    List<SnpEffInfo> getSnpEffInfo() {
        // Since we are handing back the snpeff info, 
        // it is possible the user will modify it
        this.snpEffDirty = true
        
        if(snpEffInfo != null)
            return snpEffInfo
            
        String eff = getInfo()['EFF']
        if(!eff)
            return null
            
        def genes = (eff =~ /[A-Z]*\(([^\)]*)\)/).collect { (it[1]+" ").split(/\|/)[5] }
        def txs = (eff =~ /[A-Z]*\(([^\)]*)\)/).collect { (it[1]+" ").split(/\|/)[8] }
        def effs = (eff =~ /([A-Z0-9_]*)\(([^\)]*)\)/).collect { it[1] }
        def ranks = (eff =~ /[A-Z]*\(([^\)]*)\)/).collect { (it[1]+" ").split(/\|/)[0] }
        def infos = eff.split(',')
        
        snpEffInfo = [genes,effs,ranks,infos,txs].transpose().collect { new SnpEffInfo(gene:it[0], type: it[1], impact:it[2], info: it[3], transcript:it[4]) }
    }
    
    List<Integer> getDosages() {
        if(dosages != null)
            return dosages
        
        // TODO: dosages work differently in de novo callers, 
        // they will just supply the base instead of allele 
        // Here is how it was computed in another script    
        /*
          def alleles = s.split(':')[0].split('/')
          if(alleles[0] in ['A','C','T','G']) 
              alleles[0] == alleles[1] ? 2 : 1
          else
              alleles.collect { 
                  it == '.' ? -1 : Integer.parseInt(it) 
              }.countBy { it }[1]?:0
 
         */
        dosages = genoTypes*.GT*.split('/')*.sum { it == '.' ? 0 : Integer.parseInt(it) } 
    }

    @CompileStatic 
    int sampleDosage(String sampleName) {
        return getDosages()[header.samples.indexOf(sampleName)]
    }
    
    /**
     * Returns true if this variant segregates in a way compatible with the 
     * phenotype in the given pedigree.
     * That is, if any affected sample (with a phenotype) contains the variant,
     * then all the affected samples must. If no affected samples contain 
     * the variant then this method returns true.
     * <p>
     * Note: this does not check that unaffected samples do NOT contain the 
     * variant.
     * 
     * @param p Pedigree to test 
     * @return
     */
    boolean segregatesWith(Pedigree ped) {
        // Nobody with the variant is affected OR everybody with the variant is affected
        !ped.affected.any {sampleDosage(it)} || ped.affected.every {sampleDosage(it)} 
    }
    
    /**
     * Return a list of all the pedigrees that contain this variant
     */
    Set<Pedigree> getPedigrees() {
        if(this.@pedigrees == null) {
          def dsgs = getDosages()
          this.@pedigrees = new HashSet()
          for(int i=0; i<dsgs.size(); ++i) {
              if(dsgs[i] > 0) {
                  def ped = header.findPedigreeBySampleIndex(i)   
                  if(!ped) 
                      System.err.println("WARNING: no pedigree information found for sample " + header.samples[i])
                  else
                      this.@pedigrees << ped
              }
          }
        }
        return this.@pedigrees
    }
    
    /**
     * Apply this variant to mutate the given sequence at the specified position
     */
    String applyTo(String sequence, int position) {
        switch(this.type) {
            case "SNP":
                return applySNP(sequence, position)
            case "INS":
                return applyINS(sequence, position)
            case "DEL":
                return applyDEL(sequence, position)
            default:
                throw new IllegalStateException("Unknown type of mutation: $type")
        }
    }
    
    String applySNP(String sequence, int position) {
        StringBuilder result = new StringBuilder()
        result.append(sequence.substring(0,position-1))
        result.append(alt)
        result.append(sequence.substring(position))
        return result.toString()
    }
    
    /**
     * Apply this deletion to the specified sequence at the specified position
     * <p>NOTE: this variant MUST be a deletion
     * @param sequence
     * @param position
     * @return given sequence mutated to include this deletion
     */
    String applyDEL(String sequence, int position) {
        if(this.type != "DEL")
            throw new IllegalStateException("Function can only be called for deletions")
            
        // The actual position of the deletion is the base AFTER the pos field
        sequence.substring(0,position) + sequence.substring(position + ref.size() - alt.size())
    }
    
    /**
     * Apply this deletion to the specified sequence at the specified position
     * <p>NOTE: this variant MUST be a deletion
     * @param sequence
     * @param position
     * @return given sequence mutated to include this deletion
     */
    String applyINS(String sequence, int position) {
        if(this.type != "INS")
            throw new IllegalStateException("Function can only be called for insertions")
            
        // The actual position of the insertion is the base AFTER the pos field
        // but the ALT field includes the prior base too, so go 1 backwards
        sequence.substring(0,position-1) + alt + sequence.substring(position)
    }
    
    /**
     * A debug function, displays the sequence difference and returns the 
     * new sequence
     */
    String displaySequenceDifference(String sequence, int position) {
      String updatedSequence = applyTo(sequence,position)
      if(type == "SNP") {
        [sequence,updatedSequence].each { seq ->
          println "${seq.substring(0,position-1)}[${seq.substring(position-1,position)}]${seq.substring(position)}"
        }
      }
      else 
      if(type == "DEL") {
          println sequence
          println "${sequence.substring(0,position)}${'-'*size()}${sequence.substring(position+size())}"
      }
      else 
      if(type == "INS") {
          println "${sequence.substring(0,position)}${'-'*size()}${sequence.substring(position)}"
          println updatedSequence
      }
      return updatedSequence
    }
    
    int size() {
        Math.abs(ref.size() - alt.size())
    }

    static Variant parse(String line) {
        Variant parsed = new Variant(line:line)
        parsed.parseFields(line)
        return parsed
    }
    
    String toString() {
        "$chr:$pos $ref/$alt"
    }
}
