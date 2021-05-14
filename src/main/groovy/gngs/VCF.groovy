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
 
import groovy.transform.CompileStatic;

import java.awt.event.ItemEvent;
import java.nio.MappedByteBuffer
import java.nio.channels.FileChannel.MapMode
import java.util.zip.GZIPInputStream;

import htsjdk.tribble.index.Block
import htsjdk.tribble.index.Index
import htsjdk.tribble.index.IndexFactory
import htsjdk.tribble.readers.TabixReader;

import org.codehaus.groovy.runtime.StackTraceUtils

@CompileStatic
class FormatMetaData {
    
    String id
    
    String description
    
    Class type
    
    int number = -1
}

class StopParsingVCFException extends Exception {
}

@CompileStatic
class VCFParseContext {
    final VCF vcf
    final Closure c
    final Map options 
    final Writer out 
    final Writer discardWriter 
            
    boolean flushedHeader = false
            
    public VCFParseContext(VCF vcf, Closure c, Map options, Writer out, Writer discardWriter) {
        this.vcf = vcf;
        this.c = c;
        this.options = options;
        this.out = out;
        this.discardWriter = discardWriter;
    }

    void flushHeader() {
        if(flushedHeader)
            return
                    
        if(options.updateHeader != null) {
            vcf.headerLines = (List<String>)((Closure)options.updateHeader)(vcf.headerLines)
        }
        vcf.printHeader(out)
                
        if(!(this.discardWriter.is(null))) 
            vcf.printHeader(this.discardWriter)
                
        flushedHeader = true                
    }
            
    void write(Appendable a, Variant v) {
        flushHeader()
        a.append(v.line)                
        a.append('\n')
    }
            
    void add(final Variant v) {
        if(this.out.is(null))
            this.vcf.add(v)
        else
            this.write(out,v)
    }
            
    void discard(final Variant v) {
        if(!this.discardWriter.is(null))
            this.write(discardWriter, v)
    }
    
    void close() {
        out.flush()
        out.close()
            
        if(discardWriter != null) {
            discardWriter.close()
        }
    }
}

@CompileStatic
class VCFIterator implements Iterator<Variant>, Closeable {

    Variant nextVariant = null
    
    Reader reader 
    
    VCF vcf
    
    VCFIterator(VCF vcf) {
        this.reader = Utils.createReader(vcf.fileName)
        this.vcf = vcf
        String line = vcf.readHeaders(reader)
        if(line) {
            nextVariant = Variant.parse(line)
            nextVariant.header = this.vcf
        }
    }

    @Override
    public boolean hasNext() {
        if(nextVariant == null)
            readNextVariant()
            
        if(nextVariant)
            return true

        if(!reader.is(null)) {
            reader.close()
            reader = null
        }
        return false
    }
    
    void readNextVariant() {
        String line = reader.readLine()
        if(line == null) {
            nextVariant = null
            return
        }
        Variant v = Variant.parse(line)
        v.header = this.vcf
        this.nextVariant = v
    }

    @Override
    public Variant next() {
        def result = nextVariant
        readNextVariant()
        return result
    }
    
    void close() {
        if(reader) {
            reader.close()
            reader = null
        }
    }
}


/**
 * The VCF class supports simple parsing, filtering, querying and updating of
 * VCF files.
 * <p>
 * The primary method is the {@link #parse(Closure)} method which can be
 * optionally passed a file name, stream, or if left out, will parse from
 * standard input. When a Closure is supplied, it can be used for filtering and
 * updating the VCF file as it is read. Each record is passed to the closure as
 * it is read, and if false is returned from the closure it will be discarded:
 * 
 * <pre>
 * // Let's find all the indels on chr10!
 * vcf = VCF.parse('test.vcf') { v ->
 *   v.chr == 'chr10' && (v.type == 'INS' || v.type == 'DEL')
 * }
 * </pre>
 * 
 * Various forms of information can be queried - the INFO section is parsed to a
 * map and genotype information is available, referred to as the "dosage" for each
 * sample with respect to the variant (an integer number of copies, 0,1,2). 
 * Support for SnpEff and VEP annotations is also provided via
 * dedicated properties that return a list of {@link SnpEffInfo}
 * objects describing the various annotations available. See the {@link Variant}
 * class for more information about individual variants.
 * <p>
 * Note that most information is lazily computed so the cost of parsing complex
 * fields is deferred unless they are asked for.
 * <p>
 * Updating variants requires a special procedure because related fields need to
 * be sycnhronised after modification, and it is a strong convention to add a
 * header line describing updates. An {@link #update} method provides a closure
 * based mechanism to make this straight forward:
 * <pre>
 * v.update { v.info.FOO='bar' }
 * </pre>
 * <p>
 * Limited support for "Pedigrees" is available. At the moment only grouping of
 * samples into families is supported and intra-family relationships are not
 * modeled.
 * <p>
 * The Iterable interface and "in" operator are supported, allowing all the normal
 * Groovy collection-manipulation methods.  For example, you can do
 * common iteration and membership tests:
 * <p>
 * <code>
 * for(Variant v in vcf) { 
 *   if(v in vcf2) { ... } 
 * }
 * </code>
 * or grouping, counting, etc:
 * <code>
 * vcf.countBy { it.chr } // Count of variants by chromosome
 * </code>  
 * <p>
 * If your desire is to lookup a specific region in a very large VCF file, you
 * should index the VCF file and use the {@link VCFIndex} class to query the
 * region specifically.
 * <p>
 * To stream a VCF without storing in memory, use the {@link #filter(Closure)}
 * methods which will print output, including header, directly to the output.
 *
 * @author simon.sadedin@mcri.edu.au
 */
class VCF implements Iterable<Variant> {
    
    List<String> headerLines = []
    
    String [] lastHeaderLine
    
    List<String> samples
    
    List<Variant> variants = []
    
    List<Pedigree> pedigrees
    
    List<Pedigree> samplePedigrees
    
    Map<String,List<Variant>> chrPosIndex = [:]
    
    String fileName
    
    boolean lazyLoad = false
    
    static final int FORMAT_COLUMN_INDEX = 8
    
    static final int SAMPLE_COLUMN_INDEX = 9
    
    /**
     * Cached list of VEP columns, when they exist in the VCF file.
     * Lazily populated when getVepColumns() is called.
     */
    private Map<String,String []> vepColumns = null
    
    /**
     * Lazily populated when getInfoMetaData is called
     */
    private Map<String, Map> infoMetaDatas
    
    /**
     * Lazily populated when getFormatMetaData is called
     */
    private Map<String,FormatMetaData> formatMetaData
    
    VCF() {
    }
    
    /**
     * Creates an empty VCF file based on the header of the given VCF file
     * <p>
     * The VCF is set to lazy mode which auto-loads the variants when some methods are called.
     * @param fileName
     */
    VCF(File file) {
        this(file.path)
    }
    
    /**
     * Creates an empty VCF file based on the header of the given VCF file
     * <p>
     * The VCF is set to lazy mode which auto-loads the variants when some methods are called.
     * @param fileName
     */
    VCF(String fileName) {
         this.fileName = fileName
         readHeadersOnly(fileName)
        lazyLoad = true
    }

    /**
     * Creates a VCF containing all the VCFs, and initialised with the header
     * from the first variant
     * 
     * @param variants
     */
    VCF(Iterable<Variant> variants) {
        if(variants.iterator().hasNext()) {
            Variant v = variants.first()
            if(v.header != null)
                this.headerLines.addAll(v.header.headerLines)
        }
        addAll(variants)
    }

    @CompileStatic
    private addAll(Iterable<Variant> variants) {
        for(Variant v in variants) {
            this.add(v)
        }
    }
    
    /**
     * Return a list of variants starting at the given location
     *
     * @param chr
     * @param pos
     * @return
     */
    @CompileStatic
    List variantsAt(String chr, int pos) {
        chrPosIndex[chr+":"+pos]
    }
    
    /*
     * Read only the headers to discover metadata about the VCF
     * 
     * @param fileName  path to VCF file
     */
    private void readHeadersOnly(String fileName) {
        Utils.reader(fileName) { Reader r ->
            readHeaders(r)
        }
    }

    /**
     * Reads the headers of the VCF from the given reader, and stores them in the 
     * headerLines field
     * 
     * @return  the first non-header line
     */
    @CompileStatic
    String readHeaders(Reader r) {
        String line = r.readLine();
        while(line != null) {
            if(line.startsWith('#')) {
                this.headerLines.add(line)
            }
            else
                break

            line = r.readLine()
        }
        parseLastHeaderLine()
        return line
    }
    
    static VCF parse(String fileName, Pedigrees peds, Closure c = null) {
        parse(fileName,peds?.families?.values() as List, c)
    }
    
    /**
     * Convenience method to accept string for parsing file
     */
    @CompileStatic
    static VCF parse(Map options=[:], String fileName, List<Pedigree> peds = null, Closure c = null) {
        if(fileName == "-")
          parse(System.in,peds,c)
        else {
            if(fileName.endsWith('.gz')) {
                (VCF)new File(fileName).withInputStream { vcfIs ->
                    vcfIs = new GZIPInputStream(vcfIs)
                    parse(options + [fileName:fileName], vcfIs,peds,c)
                }
            }
            else {
              parse(options, new File(fileName),peds,c)
            }
        }
    }
    
    static VCF parse(Map options=[:], String fileName, Closure c) {
        parse(options, fileName,null,c)
    }
    
    static VCF parse(File f, List<Pedigree> peds = null, Closure c = null) {
        parse([:], f, peds, c)
    }
    
    @CompileStatic
    static VCF parse(Map options,File f, List<Pedigree> peds = null, Closure c = null) {
        (VCF)Utils.reader(f.absolutePath) { r ->
            parseImpl(options+[fileName:f.path],r, false, peds, c)
        }
    }
    
    static VCF parse(Closure c = null) {
        parse("-",null,c)
    }
    
    static VCF parse(List<Pedigree> peds, Closure c = null) {
        parse("-",peds,c)
    }
    
    static VCF parse(InputStream f, List<Pedigree> peds = null, Closure c = null) {
        parse([:],f,peds,c)
    }
    
    static VCF parse(Map options, InputStream f, List<Pedigree> peds = null, Closure c = null) {
        parseImpl(options,f,false, peds,c)
    }
    
    static void filter(File f, Closure c = null) {
        filter(f,null,c)
    }
    
    static void filter(File f, List<Pedigree> peds, Closure c = null) {
        InputStream stream = null
        if(f.name.endsWith('.gz')) {
            stream = new BufferedInputStream(new GZIPInputStream(new FileInputStream(f)))
        }
        else {
            stream = new BufferedInputStream(new FileInputStream(f))
        }
        filter([fileName:f.path],stream,peds,c)
    }
    
    static void filter(String fileName, Closure c = null) {
        filter(fileName,null,c)
    }
    
    static void filter(String fileName, List<Pedigree> peds, Closure c = null) {
        if(fileName == "-")
          filter([:],(InputStream)System.in,peds,c)
        else
          filter(new File(fileName),peds,c)
    }
    
    static void filter(Closure c = null) {
        filter("-",null,c)
    }
    
    static void filter(Map options=[:], InputStream f, Closure c) {
        parseImpl(options,f,true,null, c)
    }
    
    static void filter(Map options=[:], InputStream f, List<Pedigree> peds = null, Closure c = null) {
        parseImpl(options,f,true,peds,c)
    }
    
    @CompileStatic
    void each(Closure c) {
        if(lazyLoad) {
            load { v ->
                c(v)
                false
            }
        }
        else {
            this.iterator().each(c)
        }
    }
    
    VCF load(Closure c = null) {
        VCF self = this
        Utils.createStream(this.fileName).withStream { ins ->
            VCF.parse(ins, null, c, vcf:self)
        }
        lazyLoad = false
        return this
    }
    
    /**
     * Attempts to identify which version of human genome build this VCF
     * was created from. Not applicable to non-human genomes.
     * 
     * @return  hg19,hg18,GRCh37,GRCh38,hg38 or null if one of these is not identified
     */
    @CompileStatic
    String sniffGenomeBuild() {
        String line = this.headerLines.find {it.contains('##reference')}
        
        if(!line)
            return null

        for(build in ["GRCh37","GRCh38","hg18","hg19","hg38",]) {
            if(line.contains(build))
                return build
        }
        
        if(line.contains('assembly38')) {
            return 'GRCh38'
        }
        else 
        if(line.contains('assembly37')) {
            return 'GRCh37'
        }
    }
    
    /**
     * Parse the given VCF file, returning a full VCF, optionally filtered
     * by the closure (only rows returning true from closure will be included)
     * <p>
     * Several optional parameters can be passed in the options map:
     * <li>updateHeader - a closure that will be passed the final List<String> of header lines which can be
     *                    mutated when filtering
     * <li>fileName     - set or override the automatic fileName associated to the VCF
     * <li>vcf          - associate variants to given VCF instead of creating a new, empty one
     * <li>writer       - when filtering, write to given writer instead of stdout
     * <li>outputStream - when filtering, write to given output stream instead of stdout
     *
     * @param f file to parse
     * @param peds List of pedigree objects describing sample - family relationships
     * @param c filter closure
     */
    @CompileStatic
    private static VCF parseImpl(Map options=[:], InputStream f, boolean filterMode, List<Pedigree> peds, Closure c) {
        Reader r = f.newReader()
        try {
            return parseImpl(options, r, filterMode, peds, c)
        }
        finally {
            r.close()
        }
    }
    
    @CompileStatic
    static VCF filter(Map options, String fileName, Closure c) {
        new File(fileName).withReader { r ->
            parseImpl(options, r, true, null, c)
        }
    }
 
    
    @CompileStatic
    static VCF filter(Map options=[:], Reader r, Closure c) {
        parseImpl(options, r, true, null, c)
    }
    
    @CompileStatic
    static VCF parse(Map options=[:], Reader r, boolean filterMode, Closure c) {
        parseImpl(options, r, filterMode, null, c)
    }
     
    
    @CompileStatic
    private static VCF parseImpl(Map options=[:], Reader r, boolean filterMode, List<Pedigree> peds, Closure c) {
        
        boolean ignoreHomRef = !(options.includeHomRef ?: false)
        
        VCF vcf = options.vcf ? (VCF)options.vcf : new VCF(pedigrees:peds)
        
        Writer out = createOutputWriter(options, filterMode)
        Writer discardWriter  = options.discardWriter ? Utils.writer((String)options.discardWriter) : null
        
        if(options.fileName)
            vcf.fileName = options.fileName
        
        String lastHeaderLine
        long lastPrintTimeMs = System.currentTimeMillis()
        int count = 0
        boolean flushedHeader = false
        Map chrPosIndex = vcf.chrPosIndex
        ProgressCounter progress 
        if(options.containsKey('progress')) {
            if(options['progress'])
                progress = (ProgressCounter)options['progress']
        }
        else {
            progress = new ProgressCounter(withTime:true)
        } 
        
        setParsingDelegate(c)
        
        List<String> samples = (List<String>)(options.samples != null ? options.samples : null)
        List<Integer> keepColumns = null
        
        Closure parseLastHeader =  {
            if(vcf.lastHeaderLine == null) {
                    
                vcf.parseLastHeaderLine()
                    
                // Modify the header to include only samples selected
                if(samples) {
                    List headerFields = vcf.headerLines[-1].tokenize('\t')
                    keepColumns = (0..8) + (List<Integer>)vcf.samples.findIndexValues{it in samples}.collect {it + 9}
                    vcf.headerLines[-1] = headerFields[keepColumns].join("\t")
                    vcf.parseLastHeaderLine()
                }
            }
        }
        
        VCFParseContext ctx = new VCFParseContext(vcf, c, options, out, discardWriter)
        
        try {
               
            r.eachLine { String line ->
                
                ++count
                
                if(progress != null)
                    progress.count()
               
                if(line.startsWith('#')) {
                    // Only add the header if the user didn't already provde a VCF
                    // This allows a customised header to be used.
                    if(!options.vcf)
                        vcf.headerLines.add(line)
                    return
                }
                
                parseLastHeader()
                
                if(keepColumns) {
                    line = line.tokenize('\t')[keepColumns].grep { it != null }.join("\t")
                }
                
                try {
                  Variant v = Variant.parse(line,ignoreHomRef)
                  if(v == null)
                      return
                      
                  v.header = vcf
                        
                  if((samples==null) || samples.any {v.sampleDosage(it)}) {
                      processParsedVariant(ctx, v)
                  }
                }
                catch(StopParsingVCFException stop) {
                    // rethrow
                    throw stop
                }
                catch(Exception e) {
                   throw new RuntimeException("Procesing failed at line $count :\n" + line, e)
                }
            }
        }
        catch(StopParsingVCFException stop) {
            // ignore
        }
        finally {
            if(vcf.lastHeaderLine == null)
                parseLastHeader()
        }
        
        if(filterMode) {
            ctx.flushHeader()

            ctx.close()
        }
        return vcf
    }

    /**
     * Set a delegate on the given closure that adds a <code>stop</code> method
     * to break parsing.
     */
    @CompileStatic
    private static setParsingDelegate(Closure c) {
        if(c != null) {
            c.delegate = new Object() {
                def stop() { throw new StopParsingVCFException() }

                Object getStop() {throw new StopParsingVCFException() }
            }
        }
    }

    /*
     * Creats an output writer depending on the filter mode and options
     */
    @CompileStatic
    private static Writer createOutputWriter(Map options, boolean filterMode) {
        if(filterMode) {
            if(options.writer) {
                return (Writer)options.writer
            }
            else
            if(options.outputStream) {
                return ((OutputStream)options.outputStream).newWriter()
            }
            else {
                return System.out.newWriter()
            }
        }
        return null
    }
    
    
    @CompileStatic
    static void processParsedVariant(final VCFParseContext ctx, final Variant v) {
        
        boolean include = true
        
        if(!ctx.c.is(null)) {
            include = !(((Closure)ctx.c)(v) == false)
        }
        
        if(include) {
            ctx.add(v)
        }
        else
            ctx.discard(v)
    }
    
    /**
     * Add the given variant to this VCF
     */
    @CompileStatic
    VCF add(Variant v) {
        v.header = this
        variants << v
        String key = v.chr+":"+v.pos
        List chrPosVars = chrPosIndex[key]
        if(!chrPosVars)
            chrPosIndex[key] = [v]
        else
            chrPosVars.add(v)
        return this
    }
    
    /**
     * Replace the sample ids in the VCF with the given ones
     * 
     * @param sampleIds
     */
    void replaceSamples(List<String> sampleIds) {
        this.samples = sampleIds.collect { it }
        this.lastHeaderLine = this.lastHeaderLine[0..8] + this.samples.join("\t")
        this.headerLines[-1] = this.lastHeaderLine.join("\t")
    }
    
    /**
     * Change the sample ids for this VCF
     * 
     * @param sampleIds
     */
    void renameSample(String fromId, String toId) {
        
        if(!this.samples.contains(fromId))
            throw new IllegalArgumentException("VCF does not contain sample $fromId")
        
        this.headerLines[-1] = (this.lastHeaderLine[0..<SAMPLE_COLUMN_INDEX] + this.lastHeaderLine[SAMPLE_COLUMN_INDEX..-1].collect {  s ->
            if(s == fromId)
                return toId
            else
                return s
        }).join("\t")
        parseLastHeaderLine()
    }
    
    /**
     * Perform a simplistic (emphasis on simplistic) merge between VCFs.
     * <p>
     * <li>Only works if the genotype fields are identical
     * <li>Only works if the variants are normalised (primitised) to 1 per line - unnormalised variants
     *     will produce unpredictable results
     * @param other
     * @return  A VCF containing all the variants from both VCFs with genotype information merged where
     *          the same variant is in both VCFs, and null genotype information when it is missing in one
     *          or both VCFs.
     */
    VCF merge(VCF other) {
        
        // Start by copying this VCF
        Regions allVariants = this.toBED()
        for(Region otherRegion in other.toBED()) {
            allVariants.addRegion(otherRegion)
        }
        
        VCF result = new VCF()
        
        // Add the samples from the other VCF to the last line
        String newLastLine = this.headerLines[-1] + "\t" + other.samples.collect { 
            (it in this.samples) ? it + "_2" : it
        }.join("\t")
        result.headerLines = this.headerLines[1..-2] + [ newLastLine ]
        result.parseLastHeaderLine()
        
        
        this.mergeDifferentVepHeaders(this, other, result)
        this.mergeDifferentVepHeaders(other, this, result)
       
        List<String> samples = this.samples + other.samples
        for(Region r in allVariants) {
            
            Variant v = r.extra
            
            if(result.find(v))
                continue
                
            List fields = (v.line.split("\t") as List)
            
            // If the variant is already in this VCF, then just add the sample genotype / dosage
            Variant myVariant = this.find(v)
            Variant otherVariant = other.find(v)
            
            String [] mySplit = myVariant?.line?.tokenize("\t")
            String [] otherSplit = otherVariant?.line?.tokenize("\t")
            List<String> myGtFields 
            List<String> otherGtFields 
            
            
            List commonGtFields
            List<String> allGTFields = []
            if(myVariant && !otherVariant) {
                myGtFields = mySplit[8].tokenize(':')
                commonGtFields = myGtFields
                allGTFields.addAll(myGtFields)
            }
            else
            if(!myVariant && otherVariant) {
                otherGtFields = otherSplit[8].tokenize(':')
                commonGtFields = otherGtFields
                allGTFields.addAll(otherGtFields)
            }
            else {
                // Variant in both - find the common set
                if(mySplit[8] == otherSplit[8]) { // hopefully this is most of the time!
                    myGtFields = otherGtFields = commonGtFields = mySplit[8].tokenize(':')
                    allGTFields.addAll(myGtFields)
                }
                else {
                    myGtFields = mySplit[8].tokenize(':')
                    otherGtFields = otherSplit[8].tokenize(':')
                    commonGtFields = myGtFields.intersect(otherGtFields)
                }
                allGTFields = (myGtFields + otherGtFields).unique()
            }
            
//            int numGtFields = commonGtFields.size()
            int numGtFields = allGTFields.size()
            
            List newLine = fields[0..7] + [allGTFields.join(":")]
                
            if(myVariant) {
                // The variant is in my VCF - use the genotypes from there
                def gts = mySplit[SAMPLE_COLUMN_INDEX..-1]
                newLine.addAll(gts.collect { gt ->
                    buildGtFieldValue(gt, myGtFields, allGTFields)
                })
            }
            else { 
                // Add null genotypes from our sample
                newLine.addAll(this.samples.collect { (["0/0"] + ["."] * (numGtFields-1)).join(':') })
            }
            
            if(otherVariant) {
                // The variant is in the other VCF - use the genotypes from there
                // TODO: if different number of alleles in other sample, this will not work!
                String [] otherSplit2 = otherVariant.line.split("\t")
                if(otherSplit2.size()<=SAMPLE_COLUMN_INDEX)
                    print "WTF?"
                def gts = otherSplit2[SAMPLE_COLUMN_INDEX..-1]
                newLine.addAll(gts.collect { gt ->
                    buildGtFieldValue(gt, otherGtFields, allGTFields)
                }) 
            }
            else {
                // Add null genotypes from the other samples
                newLine.addAll(other.samples.collect { (["0/0"] + ["."] * (numGtFields-1)).join(":") })
            }
            
            Variant resultVariant 
            try {
                resultVariant = Variant.parse(newLine.join("\t"))
            }
            catch(Exception e) {
                throw new IllegalStateException("Error merging variant records. Merged line = " + newLine, e)
            }
            
            result.add(resultVariant)
        }
        return result
    }
    
    @CompileStatic
    String buildGtFieldValue(String gt, List<String> myGtFields, List<String> includedFields) {
        Map<String,String> gtFields = [myGtFields,gt.tokenize(":")].transpose().collectEntries()
        includedFields.collect { String gtField ->
            String value = gtFields[gtField]
            if(!value)
                value = '.'
            return value
        }.grep { it != null }.join(":") // note: filtering by null necessary because the format field can be
                                        // longer than the gt, eg: GT:DP  ./.
    }
    
    /**
     * Check if the first vcf uses a CSQ VEP info tag while the second uses
     * an ANN info tag. If so, add the ANN info tag to the result.
     * 
     * this is an internal function for solving a specific case of merging
     * two VCFs with VEP annotations done using different VEP versions / annotation
     * info tags.
     * 
     * @param vcf1
     * @param vcf2
     * @param result VCF to add the ANN tag to
     */
    private static mergeDifferentVepHeaders(VCF vcf1, VCF vcf2, VCF result) {
        
        // If the other VCF has VEP annotations in different format,
        // add an info line for the different format to our own
        if(vcf1.getInfoMetaData("CSQ") && vcf2.getInfoMetaData("ANN")) {
            
            // Find where to add
            int insertionPoint = result.headerLines.findLastIndexOf { it.startsWith('##INFO=') }
            
            // Find line to add
            String lineToAdd = vcf2.headerLines.find { it.startsWith('##INFO=<ID=ANN,Number') }
            if(!lineToAdd)
                throw new Exception("Unable to locate INFO definition for VEP annotations (ANN) in VCF to be merged")
            
            result.headerLines.add(insertionPoint, lineToAdd)
        }
 
    }
    
    Pedigree findPedigreeBySampleIndex(int i) {
        def peds = getPedigrees()
         if(samplePedigrees == null) {
             samplePedigrees = []
             samples.eachWithIndex { String s, int index ->
                 samplePedigrees[index] =  peds.find { p -> p.individuals.any { it.id == s } }
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
     * 
     * If samples are specified, the samples parsed will be limited to the given
     * samples. The VCF header line will be modified so that any output will only 
     * include the given samples.
     */
    @CompileStatic
    void parseLastHeaderLine() {
        this.lastHeaderLine = this.headerLines[-1].split('[\t ]{1,}')
        this.samples = []
        if(this.lastHeaderLine.size()<=SAMPLE_COLUMN_INDEX) {
            this.samples = ['NA']
            if(this.lastHeaderLine.size()<=FORMAT_COLUMN_INDEX) {
                this.lastHeaderLine = this.lastHeaderLine + ['FORMAT','NA']
            }
            else {
                this.lastHeaderLine = this.lastHeaderLine + ['NA']
            }
            this.headerLines[-1] = this.lastHeaderLine.join('\t')
        }
        else
        for(int i=SAMPLE_COLUMN_INDEX; i<this.lastHeaderLine.size(); ++i) {
            this.samples.add(this.lastHeaderLine[i])
        }
    } 
    
    @CompileStatic
    List<String> getContigs() {
    
        // ##contig=<ID=chr17_gl000203_random,length=37498,assembly=hg19>
        
        headerLines.grep { String line -> line.startsWith('##contig=') }
                   .collect { String line ->
                       int idIndex = line.indexOf('ID=')
                       int commaIndex = line.indexOf(',', idIndex)
                       return line.substring(idIndex, commaIndex).tokenize('=')[1].trim()
                   }.grep { it }
    }
    
    @CompileStatic
    Map<String,Object> getInfoMetaData(String id) {
        if(this.infoMetaDatas == null) {
          this.infoMetaDatas = 
              this.headerLines.grep { String hline -> hline.startsWith("##INFO") }
                              .collect { String hline -> parseInfoMetaData(hline) }
                              .collectEntries { Map hEntry -> 
                                  [ hEntry.ID, hEntry ] 
                              }
        }
        return this.infoMetaDatas[id]
    }
    
    Map<String,FormatMetaData> getFormatMetaData() {
        if(this.formatMetaData == null) {
            this.formatMetaData = 
              this.headerLines.grep { String hline -> hline.startsWith("##FORMAT") }.collect { String hline -> parseFormatMetaDataLine(hline) }.collectEntries { 
                  FormatMetaData metaData -> [ metaData.id, metaData ] 
              }            
        }
        
        return this.formatMetaData
    }
    
    FormatMetaData parseFormatMetaDataLine(String line) {
        String attString = line[10..-2]
        
        Map attributes = attString.tokenize(',').collectEntries { String attList ->
            attList.tokenize('=')
        }
        
        Class type
        switch(attributes.Type) {
            case "String":
                type = String
                break
                
            case "Integer":
                type = Integer
                break
                
            case "Float":
                type = Float
                break
        }
        
        new FormatMetaData(type:type, description:attributes.Description[1..-2], id: attributes.ID, number: attributes.Number.isNumber()?attributes.Number.toInteger():-1)
    }
    
    /**
     * Add a header for describing an info value to be added to a 
     * VCF. 
     * 
     * @param id
     * @param desc
     * @param value A prototype value to infer the type of value for the INFO field from.
     */
    void addInfoHeader(String id, String desc, Object value) {
        int lastInfo = this.headerLines.findLastIndexOf { it.startsWith("##INFO=") }
        if(lastInfo < 0)
            lastInfo = 1
                    
        String valueType = "String"
        if(value instanceof Integer) {
            valueType = "Integer"
        }
        else
        if(value instanceof Float) {
            valueType = "Float"
        }
        else
        if(value instanceof Double) {
            valueType = "Double"
        }
                     
        this.headerLines = this.headerLines[0..<lastInfo] + 
            ["##INFO=<ID=$id,Number=1,Type=${valueType},Description=\"$desc\">"] +
            this.headerLines[lastInfo..-1]
                        
        // Clear any cached meta data so it will get reparsed
        this.infoMetaDatas = null;
    }
    
    /**
     * Return true if this VCF file contains the specified INFO tags
     * <p>
     * Note: it only checks if the tag is described in the header,
     * not whether any variant in the VCF actually has the INFO tag.
     * You still need to account that any given record may be missing the
     * tag.
     * 
     * @param id    id of INFO tag to check for
     * @return  true iff an INFO meta data line of the VCF file describes this tag
     */
    @CompileStatic
    boolean hasInfo(String id) {
        return getInfoMetaData(id) != null
    }
    
    @CompileStatic
    boolean hasVEP() {
        return hasInfo('CSQ') ||  // standard
               hasInfo('ANN') ||  // adopted by MGHA
               hasInfo('vep') // gnomAD
    }
    
    @CompileStatic
    String[] getVepColumns() {
        getVepColumns("auto")
    }
    
    /**
     * Dedicated support for dynamically returning the VEP columns present in this VCF file
     * 
     * @return
     */
    @CompileStatic
    String[] getVepColumns(String vepType) {
        parseVepColumns()
        
        if(vepType == "auto") {
            return this.vepColumns.values()[0]
        }
        else {
            return this.vepColumns[vepType]
        }
    }
    
    Map<String,String[]> parseVepColumns() {
        
        if(vepColumns != null)
            return vepColumns
            
        vepColumns = [:]
        
        for(vepType in ["CSQ","ANN","vep"]) { 
               
            Map csqInfo = getInfoMetaData(vepType)
               
            if(!csqInfo)
                continue

            def vepFormatMatch = (csqInfo.Description =~ 'Format: (.*)')
            if(!vepFormatMatch)
                continue

            String fields = vepFormatMatch[0][1]
            
            String [] vepTypeColumns = fields.split("\\|")
            
            vepColumns[vepType] = vepTypeColumns
        }
        
        if(vepColumns.size() == 0) 
            throw new IllegalArgumentException("VCF does not contain VEP annotations")
        
        return vepColumns
    }
    
    /**
     * Parse a VCF INFO meta data line and return the values as a Map
     * 
     * @param info
     * @return  map containing entries: ID, Number, Type, Description
     */
    Map<String, Object> parseInfoMetaData(String info) {
        
        String contents = (info =~ /##INFO=<(.*)>/)[0][1]
            
        // remove description
        int descriptionIndex = contents.indexOf('Description=')
        String description = contents.substring(descriptionIndex)
            
        // Get index of all the fields prior to description
        Map metaFields = 
            Variant.COMMA_SPLIT.split(contents.substring(0,descriptionIndex))
                               .collectEntries {
                                   [it.tokenize("=")[0], it.tokenize("=")[1]] 
                               }
                               
              
        // If description is in quotes, strip them off
        if(description.startsWith('"') && description.endsWith('"'))
            description = description[1..-1]
            
        metaFields.Description = description
            
        return metaFields
    }
    
    // Backwards compatibility ...
    Regions toBED() {
        toRegions()
    }
    
    @CompileStatic
    Regions toRegions() {
        Regions bed = new Regions()
        for(Variant v in this) {
          bed.addRegion(v.chr, v.pos, v.pos + v.size(), v)
        }
        return bed
    }
    
    Object asType(Class clazz) {
        if(clazz == Regions || clazz == BED ) {
            return toRegions()
        }
        else {
            return super.asType(clazz)
        }
    }
    
    /**
     * Support for the convenience syntax to check if a VCF contains a specific variant
     * with the form:
     * <pre>
     * if(variant in vcf) {
     *     ....
     * }
     * </pre>
     * Note however that dosage of the variant is NOT compared, so a variant is "in" a VCF
     * even if its allele count / zygosity is different to the variant of interest.
     * 
     * @param obj
     * @return  true iff the object is a Variant and a variant with the same allele and position 
     *          exists in this VCF
     */
    @CompileStatic
    boolean isCase(Object obj) {
        if(obj instanceof Variant) {
            Variant v = (Variant)obj
            String chrPos = v.chr + ':' + v.pos
            List<Variant> variantsAtPos = chrPosIndex[chrPos]
            
            // Note that it's "in" there if any of the ALT's match, we don't require all of them
            return variantsAtPos?.any { 
                it.ref == v.ref && it.alts.any { it == v.alt }
            }
        }
    }
    
    
    /**
     * Attempt to locate a variant having the same change as the given variant
     * inside this VCF.
     * <p>
     * This function does <b>not</b> guarantee to find the change, if it exists.
     * Currently, it will only locate the change if it is represented by an entry 
     * in the VCF starting at the same reference position. A variant will be returned
     * if <b>any</b> allele in the identified variant matches <b>any</b> allele in the
     * given variant. These operations become much more reliable (but still not
     * 100% reliable) when both VCFs from which the variants are drawn are decomposed
     * into primitives.
     * <p>
     * NOTE: this function does not require a genotype level match to return a variant. It
     * will return a variant found in the VCF regardless of genotype (eg: even hom ref)
     * 
     * @param v Variant to find
     * @return a Variant matching v, if one is found, otherwise null
     */
    @CompileStatic
    Variant find(Variant v) {
        List<Variant> candidates = chrPosIndex[v.chr+':'+v.pos]
        if(!candidates)
            return null
            
        candidates.grep { Variant vc -> vc.ref == v.ref }.find { 
            it.alleles.any { myAllele -> 
                v.alleles.any { 
                    otherAllele -> myAllele.alt == otherAllele.alt && myAllele.start == otherAllele.start  
                } 
            }
        }
    }
    
    /**
     * Calculate the number of variants for each different VEP consequence, returning
     * a map of consequence => count.
     */
    Map<String, Integer> getConsequenceCounts() {
        this.collect { v -> 
            int i=0; 
            v.alleles.collect { a -> v.getConsequence(i++) } } // get consequences for all alleles as a list
                     .flatten() // flatten the list to get a flat list of all consequences
                     .collect { it?:"unknown" } // Some nulls appear, not sure why
                     .countBy { it } // count by consequence
    }
    
    /**
     * Add new INFO header fields to this VCF
     * <p>
     * Searches for the last position of the other INFO lines and inserts the new lines at that position.
     * If there are no other INFO lines, inserts the new lines at the end of the headers
     * 
     * @param lines
     */
    void addInfoHeaders(List<String> linesToAdd) {
        // Add format field
        int lastInfoIndex = headerLines.findLastIndexOf { it.startsWith('##INFO') }
        if(lastInfoIndex < 0)
            lastInfoIndex = headerLines.size()-2
        
        this.headerLines = headerLines[0..lastInfoIndex] + 
                linesToAdd + 
                headerLines[(lastInfoIndex+1)..-1]
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
                    genes[eff.gene] = (List<Variant>)[]
                }
                genes[eff.gene] << v
              }
          }
        }
        return genes
    }
    
    /**
     * Print out a version of this VCF with only the given samples included
     */
    void filterSamples(PrintStream p, List<String> includeSamples) {
        
        p.println(headerLines[0..-2].join('\n'))
        
        // Remove the unwanted samples from the last header line
        List<String> lastHeaderFields = lastHeaderLine[0..8] + includeSamples
        p.println lastHeaderFields.join("\t")
        def indices = includeSamples.collect { s -> samples.indexOf(s) }*.plus(SAMPLE_COLUMN_INDEX)
        System.err.println "Indices of selected samples are $indices"
        for(Variant v in this) {
            // Only write the variant if it has non-zero dosage for at least one sample
            if(includeSamples.any { v.sampleDosage(it)}) {
                String [] variantFields = v.line.split("\t")
                p.println((variantFields[0..8] + variantFields[indices]).join("\t"))
            }
        }
    }
    
    /**
     * Finds a set of high quality SNPs that distinguish this VCF from the
     * other VCF provided.
     */
    List<Variant> getHighQualityHets(VCF p2) {
        return this.grep { Variant v ->
            (v.type == "SNP") && 
            (v.qual > 20) &&
            (v.dosages[0] == 1) &&
            (v.alleleBalance < 0.1) &&
            !p2.variantsAt(v.chr, v.pos)
        }        
    }
    
    /**
     * Calculate the rate at which variants are transmitted to this sample from p1,
     * given p2 is the other parent.
     * 
     * @param parent1
     * @param parent2
     * @return
     */
    double transmissionRate(VCF p1, VCF p2) {
//        def counts = (p1.chrPosIndex.keySet() - p2.chrPosIndex.keySet()).countBy { it in chrPosIndex.keySet() } /* +
//                     (p2.chrPosIndex.keySet() - p1.chrPosIndex.keySet()).countBy { it in this.chrPosIndex.keySet() } */
//            
        List<Variant> markers = p1.getHighQualityHets(p2)
        
        Map<Boolean,Integer> counts = markers.countBy { Variant v ->
            variantsAt(v.chr, v.pos) != null
        }
        
        if(!counts[true]) {
            return 0.0d
        }
        
        if(!counts[false]) {
            return Double.MAX_VALUE
        }
  
        return counts[true] / (Double)(counts[false] + counts[true])
    }
    
    /**
     * Infer proband and calculate rate of de novo variants for a VCF containing a trio.
     * 
     * @return Map in the form [proband: <sample id>, denovoRate: <fraction of denovos>, ratio: <fraction of denovos / min of other two>]
     */
    Map trioDenovoRate() {
        
        List combos = [
            [0,1,2],
            [1,0,2],
            [2,0,1]
        ]
        
        int i=0
        List<Map> rates = combos.collect { List combo ->
            [ 
                proband: samples[combo[0]],
                denovoRate: denovoRate(*combo),
                index: i++
            ]
        }.sort {
            it.denovoRate
        }
        
        Map result = [
            proband: rates[0].proband,
            denovoRate: rates[0].denovoRate,
            ratio: rates[0].denovoRate / rates[1].denovoRate,
            trios: rates
        ]
        
        return result
    }

    /**
     * Calculate the fraction of denovo variants if this VCF contains a trio with
     * the indexes of the samples as given
     * 
     * @param proband   index in the samples array of the proband
     * @param parent1   index in the samples array of first parent
     * @param parent2   index in the samples array of second parent
     * @return
     */
    @CompileStatic
    double denovoRate(int proband, int parent1, int parent2) {
        Map<Boolean,Integer> counts = this.grep { Variant v ->
            (v.type == "SNP") && (v.dosages[proband] == 1) && (v.qual > 20) 
        }
        .countBy { Object o ->
            Variant v = (Variant)o
            List ds = v.dosages
            
            // true of variant IS de novo
            (ds[proband]>0) && ((ds[parent1] == 0)  && (ds[parent2] == 0))
        }
        
        // true => de novo
        // false => inherited

        if(!counts[false]) {
            return 1.0d // no inherited => rate = 100% denovo
        }
        
        if(!counts[true]) {
            return 1.0d // no de novo => rate = 0% de novo
        }
  
        return counts[true] / (double)(counts[false] + counts[true])
    }

    double denovoRate(VCF m, VCF d) {
        Map<Boolean,Integer> counts = this.grep { Variant v ->
            (v.type == "SNP") && (v.dosages[0] == 1) && (v.qual > 20) 
        }.countBy { Object o -> Variant v = (Variant)(o); (v in m) || (v in d) }
        
        if(!counts[true]) {
            return 1.0d
        }
        
        if(!counts[false]) {
            return 0.0d
        }
  
        return counts[false] / (double)(counts[false] + counts[true])
    }
    
    /**
     * Return the sex of a sample, estimated from the heterozygosity of its variants.
     * <p>
     * See {@link #guessSex(int)} for details of the implementation.
     * 
     * @param sampleId
     * @return  estimated Sex of Sample
     */
    Sex guessSex(String sampleId, int sampleSize=500) {
        if(!this.samples.contains(sampleId))
            throw new IllegalArgumentException("Sample $sampleId is not contained in this VCF")
            
        guessSex(this.samples.indexOf(sampleId), sampleSize)
    }
    
    
    /**
     * Attempt to guess the sex of a human VCF file by sampling high quality 
     * variants on the X chromosome.
     * <p>
     * With various checks and constraints, the main test is: 
     *   <li> if the number of homozygotes > 2 * number of hets,
     *        then sex = FEMALE
     *   <li> else sex = MALE
     *   <li> if various checks fail (eg: less than 100 variants available),
     *        then this method may return {@link Sex.UNKNOWN}.
     *   
     * <em>Note: requires an indexed VCF file</em>
     * 
     * @param vcf
     * @return  estimated Sex of sample
     */
    Sex guessSex(int sampleIndex = 0, int sampleSize=500) {
        
        final int minDepth = 13
        
        final int minFinalVariantsForSexEstimation = 50
        
        // open the VCF and sample 100 variants from the X and Y chromosomes
        VCFIndex index = new VCFIndex(fileName)
        
        String chr = 'chrX'
        if(new File(fileName + '.tbi').exists()) {
            TabixReader r = new TabixReader(fileName, fileName + '.tbi')
            if('X' in r.chromosomes) {
                chr = 'X'
            }
        }
        
        // This region is stable in both GRCh38 and GRCh37 for sampling 
        final Region variantSamplingRegion = new Region(chr,5000000,155270560)
            
//        else {
//            System.err.println "WARNING: No contig index - assuming chrX for sex chromosome"
//        }
        
        int sexEstimationVariantCount = sampleSize
        
        List<Variant> xVariants = new ArrayList(sexEstimationVariantCount+1)
        
        index.query(variantSamplingRegion) { Variant v ->
            if(xVariants.size() >= sexEstimationVariantCount)
                return false
                
            // Consider only high quality variants that do not have a highly skewed 
            // allele balance    
            int dp = v.genoTypes[sampleIndex].DP ?: 20
            if(dp  >= minDepth && v.qual > 20 && v.alleleBalance<0.1)
                xVariants << v
        }
        
        if(xVariants.size() < minFinalVariantsForSexEstimation)
            throw new IllegalStateException("Too few variants available to use for estimation")
        
        Map<Integer, Integer> dosages = xVariants.countBy { it.getDosages(0)[sampleIndex] }
        
        // No heterozygotes
        if(dosages[1] == null) {
            return Sex.MALE
        }
        
        // No homozygotes 
        if(dosages[2] == null) {
            if(xVariants.size()>100) // make sure we sampled enough
                return Sex.FEMALE
            else
                return Sex.UNKNOWN
        }
         
        // Excess of homozygous?
        if(dosages[2] > dosages[1] * 2)
            return Sex.MALE
        else
            return Sex.FEMALE
    }
    
    void printHeader(Appendable w) {
        w.append(headerLines.join('\n'))
        w.append('\n')
    }
  
    void printHeader() {
        printHeader(System.out)
    }
    
    void print() {
        this.print(System.out)
    }
    
    @CompileStatic
    void print(PrintWriter p) {
        print((Appendable)p)
    }
    
    @CompileStatic
    void print(Appendable p) {
        printHeader(p)
        for(Variant v in variants) {
            p.append(v.line)
            p.append('\n')
        }
    }

    /**
     * Return an iterator that iterates through the variants in this
     * VCF.
     * 
     * <p><b>Note</b>: if the VCF is created lazily (file name constructor) then the 
     * iterator must be closed if iteration does not complete.</p>
     */
    @Override
    @CompileStatic
    public Iterator<Variant> iterator() {
        if(this.lazyLoad) {
            return new VCFIterator(this)
        }
        else
            return this.variants.iterator()
    }
    
    int getSize() {
        this.variants.size()
    }
    
    @CompileStatic
    List<Map<String,Object>> toListMap() {
     (List<Map<String,Object>>)[
            this.toRegions().toListMap(), 
            this*.info, 
            this*.genoTypes*.getAt(0),
            this*.ref,
            this*.alt,
        ].transpose().collect { Object infosObj -> List infos = (List)infosObj;  ((Map)infos[0]) +  [Ref:infos[3], Alt: infos[4]] + ((Map)infos[1]) + ((Map)infos[2])}
    }
    
    private VCFIndex index
    
    VCFIndex getIndex() {
        if(!index)
            index = new VCFIndex(fileName)
        return index
    }
    
    String toString() {
        "VCF ${this.fileName != null ? ('file ' + this.fileName) : '' } for ${samples?.join(',')}"
    }
}