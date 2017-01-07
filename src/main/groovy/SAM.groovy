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

import java.lang.management.ManagementFactory;
import java.lang.management.OperatingSystemMXBean;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import groovy.stream.Stream;

//import com.sun.management.UnixOperatingSystemMXBean;

import groovy.transform.CompileStatic;
import htsjdk.samtools.BAMRecord;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory
import htsjdk.samtools.SAMFormatException;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

@CompileStatic 
class SAMRecordPair {
    SAMRecord r1
    
    SAMRecord r2
    
    boolean hasBothReads() {
        r2 != null && r1 != null
    }
    
    boolean isChimeric() {
        return r1 != null && r2 != null && (r1.referenceIndex != r2.referenceIndex)
    }
}

/*
 * An alternate alignment for a read. 
 * eg: constructed from the XA tag as output by bwa aln.
 * 
 * @author Simon
 */
class XA {
    String chr
    int pos
    String cigar
    int nm

    public String toString() {
        "[$chr:$pos $cigar NM=$nm]"
    }
}

class SAMRecordCategory {
    @CompileStatic
    static List<XA> getAlternateAlignments(BAMRecord record) {
        String xas = record.getAttribute("XA")
        if(!xas)
            return null
        return xas.split(';').collect { String xa ->
            String [] parts = xa.split(',')
            new XA( chr: parts[0], pos: parts[1].substring(1) as Integer, cigar:parts[2],  nm : parts[-1] as Integer)
        }
    }
    
    @CompileStatic
    static Object asType(SAMRecord r, Class clazz) {
        new Region(r.referenceName, Math.min(r.alignmentStart, r.alignmentEnd)..Math.max(r.alignmentStart, r.alignmentEnd))
    }
    
    @CompileStatic
    static Object asType(BAMRecord r, Class clazz) {
        new Region(r.referenceName, Math.min(r.alignmentStart, r.alignmentEnd)..Math.max(r.alignmentStart, r.alignmentEnd))
    } 
    
    @CompileStatic
    static Object toRegion(BAMRecord r) {
        new Region(r.referenceName, Math.min(r.alignmentStart, r.alignmentEnd)..Math.max(r.alignmentStart, r.alignmentEnd))
    }  
}

@CompileStatic
class ReadWindow {
    
    int pos
    
    TreeMap<Integer, List<SAMRecord>> window = new TreeMap()
}

/**
 * Adds various Groovy idioms and convenience features to the 
 * Picard SAMFileReader.
 * <p>
 * There are three major classes of functionality supported:
 * 
 * <li>Iterating through and filtering reads in various ways
 * <li>Generating pileups and calculating read depth / coverage
 * <li>Accessing meta data (read groups, sample information, etc)
 * 
 * For simple looping, the {@link #eachRead(Closure)} static method
 * can be used without creating a SAM object at all:
 * <pre> SAM.eachRead { SAMRecord r -> println r.readName } </pre>
 * (this will read from standard input). More sophisticated use requires the construction of a SAM object, 
 * which allows, for example, iteration of read pairs:
 * <pre> new SAM("test.bam").eachPair { r1, r2 -> assert r1.readName == r2.readName } </pre>
 * <p>
 * A region can be optionally passed to iterate over:
 * 
 * <pre> new SAM("test.bam").eachPair("chr1",1000,2000) { r1, r2 -> assert r1.readName == r2.readName }</pre>
 * 
 * Filtering a BAM file to create file containing a subset of reads is supported explicitly:
 * 
 * <pre>new SAM("test.bam").filter("out.bam") { it.mappingQuality > 30 }</pre>
 * 
 * Generating pileups is also straightforward:
 * <pre> new SAM("test.bam").pileup("chr1",1000,2000) { p -> 
 *     println "There are ${p.countOf('A')} A bases at position chr1:$p.position" 
 * }</pre>
 * See the {@link Pileup} class for more information about operations on pileups.
 * <p><br>
 * <i>Notes:</i>
 * <li>Most operations filter out reads with mapping quality 0 by default. To 
 * avoid this, set the minMappingQuality property on the SAM object.
 * <li>All operations with this class require the BAM file to be indexed.
 * <li>SAM and BAM files are treated transparently the same.
 * <p>
 * @author simon.sadedin@mcri.edu.au
 */
class SAM {

    SAMFileReader samFileReader;

    File samFile

    File indexFile
    
    /**
     * Only used when created from stream
     */
    InputStream samStream

    int minMappingQuality = 1

    static boolean progress = false

    static {
        // Picard has started doing some really verbose logging - turn it off by default
        // and let the user turn it back on with this property if they want it
        if("true" != System.properties.picardLogging)
            htsjdk.samtools.util.Log.setGlobalLogLevel(htsjdk.samtools.util.Log.LogLevel.WARNING)
    }
    
    SAM(InputStream ips) {
        SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT)
        samStream = ips
        newReader()
    }

    SAM(String fileName) {
        this(new File(fileName))
    }

    SAM(File file) {
        this.samFile = file
        if(!this.samFile.exists())
            throw new FileNotFoundException("BAM/CRAM file could not be opened at ${samFile.absolutePath}")

        String indexExt;
        String fileExt;
        if(file.name.endsWith(".bam") || file.name.endsWith(".cram")) {
            indexExt = ".bai"
            fileExt = ".bam"
        }
        else
            throw new IllegalArgumentException("This class only supports BAM or CRAM files. File type $file is not supported")

        this.indexFile = new File(samFile.absolutePath + indexExt)

        if(!indexFile.exists())
            indexFile = new File(samFile.absolutePath.replaceAll("$fileExt\$",'')+indexExt)

        if(!indexFile.exists())
            throw new FileNotFoundException("Please ensure your BAM / SAM / CRAM file is indexed. File $indexFile could not be found.")

        SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT)
        this.samFileReader = newReader()

    }

    SAMFileReader newReader() {
        if(samStream == null)
            new SAMFileReader(samFile, indexFile, false)
        else
            new SAMFileReader(samStream)
    }

    /**
     * Return a new SAMFileWriter configured with the same settings as this 
     * SAM. It is the caller's responsibility to close the writer.
     * 
     * @param outputFileName    Name of file to write to
     * @return SAMFileWriter
     */
    SAMFileWriter newWriter(String outputFileName) {
        SAMFileWriterFactory f = new SAMFileWriterFactory()
        SAMFileHeader header = this.samFileReader.fileHeader
        SAMFileWriter w = f.makeBAMWriter(header, false, new File(outputFileName))
        return w
    }

    /**
     * Return a new SAMFileWriter configured with the same settings as this 
     * SAM. It is the caller's responsibility to close the writer.
     * In this version of the method the files are assumed to be 
     * <em>pre-sorted</em>*.
     * 
     * @param outputFileName    Name of file to write to
     * @return SAMFileWriter
     */
    def withWriter(String outputFileName, Closure c) {
        withWriter(outputFileName, true, c)
    }
    
    def withWriter(String outputFileName, boolean sorted, Closure c) {
        SAMFileWriterFactory f = new SAMFileWriterFactory()
        SAMFileHeader header = this.samFileReader.fileHeader
        SAMFileWriter w = f.makeBAMWriter(header, sorted, new File(outputFileName))
        try {
            return c(w)
        }
        finally {
            w.close()
        }
    }
    
    def withOrderedPairWriter(String outputFileName, boolean sorted, Closure c) {
        SAMFileWriterFactory f = new SAMFileWriterFactory()
        SAMFileHeader header = this.samFileReader.fileHeader
        SAMFileWriter w = f.makeBAMWriter(header, sorted, new File(outputFileName))
        OrderedPairWriter opw = new OrderedPairWriter(w)
        try {
            return c(opw)
        }
        finally {
            w.close()
        }
    }

    /**
     * Iterate over each record in the same file in the order they are in the file
     * @param c
     */
    void eachRecord(Closure c) {
        use(SAMRecordCategory) {
            
            SAMFileReader reader = this.newReader()
            try {
                SAMRecordIterator i = reader.iterator()
                try {
                    while(i.hasNext()) {
                        c(i.next())
                    }
                }
                finally {
                    i.close()
                }
            }
            finally {
                reader.close()
            }
        }
    }
    
    /**
     * Call the given closure for every pair of reads in a BAM file containing paired end reads.
     * <p>
     * Note: the algorithm works by keeping a running buffer of reads, and iterating through
     * the reads in order until each single read finds its mate. This means that 
     * reads having no mate accumulate in the buffer without ever being removed. Thus 
     * a large BAM file containing millions of unpaired reads could cause this method to use
     * substantial ammounts of memory.
     * 
     * @param c        Closure to call    
     */
    @CompileStatic
    void eachPair(Closure c) {
        eachPair([:],c)
    }

    /**
     * Call the given closure for every pair of reads in a BAM file containing paired end reads.
     * <p>
     * Note: the algorithm works by keeping a running buffer of reads, and iterating through
     * the reads in order until each single read finds its mate. This means that 
     * reads having no mate accumulate in the buffer without ever being removed. Thus 
     * a large BAM file containing millions of unpaired reads could cause this method to use
     * substantial ammounts of memory.
     * 
     * @param c        Closure to call
     */
    @CompileStatic
    void eachPair(Map options, Closure c) {
        SAMRecordIterator iter = samFileReader.iterator()
        eachPair(options, iter,c)
    }
    
    boolean verbose = false

    /**
     * Call the given closure for every pair of reads in a BAM file containing paired end reads.
     * <p>
     * Note: the algorithm works by keeping a running buffer of reads, and iterating through
     * the reads in order until each single read finds its mate. This means that 
     * reads having no mate accumulate in the buffer without ever being removed. Thus 
     * a large BAM file containing millions of unpaired reads could cause this method to use
     * substantial ammounts of memory.
     * 
     * @param iter    iterator to consumer reads from
     * @param c        Closure to call
     */
    @CompileStatic
    void eachPair(Map options, SAMRecordIterator iter, Closure c) {
        SAMFileReader pairReader = newReader()
        SAMFileReader randomLookupReader = newReader()
        Map<String,Integer> buffer = new HashMap()
        List<SAMRecordPair> writeSpool = new LinkedList()
        SortedMap<Integer,List<SAMRecord>> readPairs = new TreeMap<Integer,List<SAMRecord>>()
        Set<String> preprocessedReads = new HashSet()
        int pairs = 0
        
        boolean includeUnmapped = ((boolean)options?.includeUnmapped)
        
        // The spool holds a buffer that tries to pair up reads before outputting them
        // The right size for it depends on the coverage depth and the separation between
        // paired reads. It should be able to hold enough reads that most of the time
        // a read's pair is encountered before the buffer overflows. When it overflows
        // the loop below starts doing random queries (slow) to resolve the mates for reads
        // so that the pair can be output together. 
        int maxSpoolSize = 4000
        int forcedQueries = 0
        String currentChr = null
        try {
            ProgressCounter progress = new ProgressCounter(extra:{
                "Chromosome $currentChr, Spool Size=${writeSpool.size()}, Buffer Size=${buffer.size()}, ForcedQueries=$forcedQueries, Pairs=${pairs}"
            }, withRate: true)
            
            try {
                while(iter.hasNext()) {
                    SAMRecord record = (SAMRecord)iter.next();
                    
                    String readName = record.readName
                    
                    if(currentChr != null && currentChr != record.referenceName) {
                        for(SAMRecordPair readPair in writeSpool) {
                            if(readPair.hasBothReads()) {
                                c(readPair.r1, readPair.r2)
                                continue
                            }
                            
                            SAMRecord mate = this.queryMate(randomLookupReader, readPair.r1)
                           
                            if(mate != null) {
                                c(readPair.r1, mate)
                            }
                        }
                        
                        if(verbose)
                            println "XX: ${writeSpool.size()} reads ignored because mate never encountered"
                        writeSpool.clear()
                        buffer.clear()
                    }
                        
                    progress.count()
                    
//                    println "Proper pair = " + record.getProperPairFlag()
//                    println "Read paired = " + record.readPairedFlag
                        
                    // Unmapped reads are not proper pairs, so for those we do not apply the
                    // check about proper pairs. For 
                    if(includeUnmapped && (record.readUnmappedFlag || record.mateUnmappedFlag)) {
                        if(!record.getReadPairedFlag())
                            continue
                    }
                    else
                    if(!record.getReadPairedFlag() || record.getReadUnmappedFlag() || !record.getProperPairFlag())
                        continue
                        
                    if((record.referenceIndex != record.mateReferenceIndex) || record.isSecondaryOrSupplementary()) {
                        if(verbose)
                            println "Chimeric / non-primary alignment $record.readName :" + record.referenceName + ":" + record.alignmentStart
                        continue
                    }
                            
                    if(readName in preprocessedReads) {
//                        println "XX: Already processed $r1.readName"
                        preprocessedReads.remove(readName)
                        continue
                    }
                    
                    currentChr = record.referenceName
                    
                    Integer pairIndex = buffer[readName]
                    if(pairIndex == null) {
                        buffer[readName] = writeSpool.size()
                        if(verbose)
                            println "XX: Buffer read $readName $record.referenceName:$record.alignmentStart (spool position ${writeSpool.size()})"
                        writeSpool << new SAMRecordPair(r1:record)
                        continue
                    }
                    
                    if(pairIndex >= writeSpool.size())
                        pairIndex = writeSpool.size()-1
                        
                    while(writeSpool.get(pairIndex).r1.readName != readName)
                        --pairIndex
                        
                    if(verbose)
                        println "XX: Pair found for $readName at $pairIndex ($record.referenceName:$record.alignmentStart)"
                    writeSpool[pairIndex].r2 = record
                    
                    if(writeSpool.size()>maxSpoolSize) {
                        SAMRecord mate = this.queryMate(randomLookupReader,writeSpool[0].r1)
//                        assert writeSpool[0].r2 != null
                        if(mate != null) {
                            if(verbose)
                                println "XX: Pair found for ${writeSpool[0].r1.readName} (${writeSpool[0].r1.alignmentStart}, ${mate.alignmentStart}) by random lookup "
                            writeSpool[0].r2 = mate
                            preprocessedReads.add(mate.readName)
                        }
                        else {
                            SAMRecord mateless = writeSpool.remove(0).r1
                            buffer.remove(mateless.readName)
                            if(verbose)
                                println "XX: Read $mateless.readName has no mates :-("
                        }
                        ++forcedQueries
                    }
                    else
                    if(pairIndex > 0) {
                        continue
                    }
   
                    while(writeSpool[0]?.hasBothReads()) {
                        SAMRecordPair pair = writeSpool.remove(0)
                        String writeReadName = pair.r1.readName
                        if(verbose)
                            println "XX: write pair $writeReadName ($pair.r1.referenceName:$pair.r1.alignmentStart, $pair.r2.referenceName:$pair.r2.alignmentStart)"
                        c(pair.r1, pair.r2)
                        ++pairs
                        buffer.remove(writeReadName)
                    }
                }

                List<SAMRecordPair> spoolResidue = new LinkedList()
                for(SAMRecordPair pair in writeSpool) {
                    if(pair.hasBothReads() && pair.r1.referenceName == currentChr) {
                        c(pair.r1, pair.r2)
                        buffer.remove(pair.r1.readName)
                        ++pairs                        
                    }
                    else {
                        spoolResidue << pair
                    }
                }
                writeSpool = spoolResidue
            }
            finally {
                pairReader.close()
                randomLookupReader.close()
                progress.end()
            }
        }
        finally {
            iter.close()
        }
    }
    
    @CompileStatic
    SAMRecord queryMate(SAMFileReader r, SAMRecord r1) {
        try {
            return r.queryMate(r1)
        }
        catch(SAMFormatException sfe) {
            // ignore
            // happens due to BWA secondary alignments
        }
        return null
    }

    void eachPair(Region r, Closure c) {
        this.eachPair(r.chr, r.from, r.to, c)
    }

    void eachPair(String chr, Closure c) {
        this.eachPair(chr, 0, 0, c) // Note Picard convention: 0 represents start and end of reference sequence
    }

    @CompileStatic
    void eachPair(String chr, int start, int end, Closure c) {
        SAMRecordIterator<SAMRecord> iter
        
        SAMFileReader reader = this.newReader()
        if(chr)
            iter = reader.query(chr, start,end,false)
        else
            iter = reader.iterator()

        try {
            try {
                eachPair(null, iter,c)
            }
            finally {
                iter.close()
            }
        }
        finally {
            reader.close()
        }
    }

    /**
     * Count the number of read pairs in the given region
     * @param r region to count pairs over
     * @return  number of read pairs that overlap the region
     */
    int countPairs(Region r) {
        countPairs(r.chr, r.from, r.to)
    }

    @CompileStatic
    int countPairs(String chr, int start, int end) {
        int c = 0;
        eachPair(chr,start,end) { Object r1, Object r2 ->
            ++c
        }
        return c
    }

    void eachRecord(int threads, Closure c) {
        List<SAMSequenceRecord> sequences = samFileReader.fileHeader.sequenceDictionary.sequences
        ExecutorService executor = Executors.newFixedThreadPool(threads)
        sequences.each { seq ->
            executor.execute {
                def reader = new SAMFileReader(samFile, indexFile, false)
                try {
                    use(SAMRecordCategory) {
                        Iterator i = reader.query(seq.sequenceName, 1, seq.sequenceLength, false)
                        while(i.hasNext()) {
                            c(i.next())
                        }
                    }
                }
                finally {
                    reader.close()
                }
            }
        }
        executor.shutdown()
    }

    List<SAMReadGroupRecord> getReadGroups() {
        samFileReader.getFileHeader().getReadGroups()
    }

    List<String> getSamples() {
        samFileReader.getFileHeader().getReadGroups()*.sample
    }

    void filter(Closure c) {
        filter("/dev/stdout",c)
    }

    @CompileStatic
    void filter(String outputFile, Closure c) {

        SAMFileReader reader = new SAMFileReader(samFile, indexFile, false)

        SAMFileWriterFactory f = new SAMFileWriterFactory()
        SAMFileHeader header = reader.fileHeader
        SAMFileWriter w = f.makeBAMWriter(header, false, new File(outputFile))
        SAMRecordIterator i = reader.iterator()
        int count = 0
        long lastPrintMs = System.currentTimeMillis()
        try {
            while(i.hasNext()) {
                SAMRecord r = (SAMRecord)i.next()
                if(c(r) == true) {
                    w.addAlignment(r)
                }
                if(count % 1000 == 0) {
                    if(System.currentTimeMillis() - lastPrintMs > 15000) {
                        System.err.println "${new Date()} Processed $count records"
                        lastPrintMs = System.currentTimeMillis()
                    }
                }
            }
        }
        finally {
            if(i != null)
                i.close()

            if(w)
                w.close()
        }
    }

    /**
     * Return the number of mapped reads overlapping the given position
     * 
     * @param chr   the sequence name / chromosome to query
     * @param pos   the chromosomal position to query
     * 
     * @return  the number of reads overlapping the position in the file
     */    
    @CompileStatic
    int coverage(String chr, int pos, Closure c=null) {
        return coverage(this.samFileReader, chr, pos, -1, c, minMappingQuality)
    }

    @CompileStatic
    int coverage(String chr, int pos, int end, Closure c=null) {
        return coverage(this.samFileReader, chr, pos, end, c, minMappingQuality)
    }

    @CompileStatic
    float meanCoverage(String chr, int pos, int end) {
        int total = 0
        this.pileup(chr, pos, end) { PileupIterator.Pileup p ->
            total += p.alignments.size()
        }
        return ((float)total)/ (end - pos + 1)
    }

    /**
     * Create a DescriptiveStatistics object
     */
    @CompileStatic
    CoverageStats coverageStatistics(String chr, int pos, int end) {
        CoverageStats stats = new CoverageStats(10000)
        int total = 0
        ProgressCounter progress = null
        if(this.progress)
            progress = new ProgressCounter(withTime:true, withRate:true)

        this.pileup(chr, pos, end) { PileupIterator.Pileup p ->
            stats.addValue(p.alignments.size())
            if(this.progress)
                progress.count()
        }
        stats
    }

    @CompileStatic
    CoverageStats coverageStatistics(Regions regions) {
        CoverageStats stats = new CoverageStats(10000)
        int total = 0
        ProgressCounter progress = null
        if(this.progress)
            progress = new ProgressCounter(withTime:true, withRate:true)

        // Flatten the regions down in case they overlap
        Regions flattenedRegions = regions.reduce()
        for(Region region in flattenedRegions) {
            this.pileup(region.chr, region.from, region.to) { PileupIterator.Pileup p ->
                stats.addValue(p.alignments.size())
                if(this.progress)
                    progress.count()
            }
        }
        stats
    }

    Regions asType(Class clazz) {
        if(clazz == Regions) {
            return toRegions()
        }
    }

    Regions toRegions(Region overlapping) {
        toRegions(overlapping.chr, overlapping.from, overlapping.to)
    }

    /**
     * Return reads from this SAM file that overlap the specified region
     * as a Regions object - that is, as a set of genomic intervals.
     * <p>
     * Reads that are missing a start or end alignment position are 
     * omitted.
     * <p>
     * If called without passing start or end, the start / end
     * are interpreted as the beginning / end of the chromosome
     * / reference sequence respectively.
     * 
     * @return Regions object containing intervals representing all
     *            the reads in this SAM object
     */
    @CompileStatic
    Regions toRegions(String chr=null, int start=0, int end=0) {
        Regions regions = new Regions()

        SAMRecordIterator i = null
        if(chr == null)
            i = samFileReader.iterator()
        else {
            i = samFileReader.query(chr, start, end, false)
        }

        ProgressCounter progress = new ProgressCounter()
        try {
            while(i.hasNext()) {
                SAMRecord r = (SAMRecord)i.next()

                // Ignore records with no start or no end position
                if(r.alignmentStart <= 0 || r.alignmentEnd <= 0)
                    continue

                Region region = new Region(r.getReferenceName(),
                                Math.min(r.alignmentStart, r.alignmentEnd)..Math.max(r.alignmentStart, r.alignmentEnd))

                region.range = new GRange((int)region.range.from,(int)region.range.to,region)
                region.setProperty('read',r)
                regions.addRegion(region)
                progress.count()
            }
        }
        finally {
            i.close()
        }
        return regions
    }
    
    def stream(Closure c) {
        SAMRecordIterator i = null
        SAMFileReader r = newReader()
        i = r.iterator()
        use(SAMRecordCategory) {
            Stream s = Stream.from(i)
            c.delegate = s;
            try {
                c(s)
            }
            finally {
                r.close()
            }
        }
        
    }
    
    Stream getRegions() {
        getRegions(null, 0, 0)
    }
    
    Stream getRegions(String chr, int start, int end) {

        ProgressCounter progress = new ProgressCounter()
        SAMRecordIterator i = null
        if(chr == null)
            i = samFileReader.iterator()
        else {
            i = samFileReader.query(chr, start, end, false)
        }
        
        Stream.from(i).filter { SAMRecord r ->
            if(!i.hasNext())
                i.close()
            progress.count()
            r.alignmentStart <= 0 || r.alignmentEnd <= 0
        }.map { SAMRecord r ->
            Region region = new Region(r.getReferenceName(),
                            Math.min(r.alignmentStart, r.alignmentEnd)..Math.max(r.alignmentStart, r.alignmentEnd))

            region.range = new GRange((int)region.range.from,(int)region.range.to,region)
            region.setProperty('read',r)
            region
        }
    }

    /**
     * Return read pairs from this SAM file that overlap the specified region
     * as a Regions object - that is, as a set of genomic intervals. Each 
     * region returned by this method spans from the start of the 5' read
     * to the end of the 3' read.
     * <p>
     * Reads that are missing a start or end alignment position are 
     * omitted. Unpaired reads are also omitted.
     * <p>
     * If called without passing start or end, the start / end
     * are interpreted as the beginning / end of the chromosome
     * / reference sequence respectively.
     * 
     * @return Regions object containing intervals representing all
     *            the reads in this SAM object
     */
    @CompileStatic
    Regions toPairRegions(String chr=null, int start=0, int end=0, int maxSize=0) {
        Regions regions = new Regions()

        ProgressCounter progress = new ProgressCounter()
        this.eachPair(chr,start,end) { SAMRecord r1, SAMRecord r2 ->
            int [] boundaries = Utils.array(r1.alignmentStart, r1.alignmentEnd, r2.alignmentStart, r2.alignmentEnd)
            if(boundaries.any { it == 0}) {
                return
            }

            if(r1.referenceName != r2.referenceName)
                return

            Region region = new Region(r1.referenceName, Utils.min(boundaries)..Utils.max(boundaries))
            if(maxSize>0 && region.size()>maxSize) {
                return
            }

            region.range = new GRange((int)region.range.from,(int)region.range.to,region)
            region.setProperty('r1',r1)
            region.setProperty('r2',r2)
            regions.addRegion(region)
            progress.count()
        }

        return regions
    }


    /**
     * Return the number of mapped reads overlapping the given position
     * 
     * @param r     the SAMFileReader (SAM / BAM file) containing reads
     * @param chr   the sequence name / chromosome to query
     * @param pos   the chromosomal position to query
     * 
     * @return  the number of reads overlapping the position in the file
     */
    @CompileStatic
    static int coverage(SAMFileReader r, String chr, int pos, int end=-1, Closure c=null, int minMappingQuality=1) {
        if(end == -1)
            end = pos

        ProgressCounter counter = new ProgressCounter()
        SAMRecordIterator i = r.query(chr, pos, end, false);
        try {
            int count = 0;
            while(i.hasNext()) {
                SAMRecord rec = (SAMRecord)i.next();
                if(progress)
                    counter.count()
                if(rec.getMappingQuality() < minMappingQuality) {
                    continue;
                }
                if(c == null || c(rec))
                    ++count;
            }
            return count;
        }
        finally {
            i.close();
        }
    }

    PileupIterator.Pileup pileup(String chr, int pos) {
        return pileup(chr,pos,pos).next()
    }

    @CompileStatic
    void pileup(String chr, int start, int end, Closure c) {
        PileupIterator i = pileup(chr,start,end)
        try {
            while(i.hasNext()) {
                c(i.next())
            }
        }
        finally {
            i.close()
        }
    }

    @CompileStatic
    PileupIterator pileup(String chr, int start, int end) {
        PileupIterator p = new PileupIterator(new SAMFileReader(samFile, indexFile, false), chr,start,end);
        p.setMinMappingQuality(this.minMappingQuality)
        return p
    }
    
    
    /**
     * Return a set of read counts indicating the counts of reads in that overlap the target region.
     * <p>
     * 
     * @todo    this is pretty inefficient compared to how it could be. It makes 2 passes through the data
     *          where it could do 1.
     * @param targets
     * @return Map with keys on_target: reads overlapping target region, and bp_on_target: number of base pairs
     *         sequenced that overlap the target region, total: total count of reads.
     */
    Map<String,Long> countOnTarget(Regions targets) {
        
        long total = 0
        def on_target = this.stream {
          map { it as Region }.filter { ++total; targets.overlaps(it) }.count { 1 }
        }
        

        long totalBp = 0
        this.stream {
          map { it as Region }.map { targets.intersect(it)[0] }.filter { it != null }.each { totalBp += it.size() }
        }
        
        [
          total: total,
          bp_on_target: totalBp,
          on_target: on_target
        ]
    }
    
    Map<String, Integer> getContigs() {
        SAMFileReader reader = newReader()
        try {
            List<SAMSequenceRecord> sequences = reader.getFileHeader().getSequenceDictionary().getSequences()
            Map<String,Integer> contigs = sequences.collectEntries { seq ->
                [seq.sequenceName, seq.sequenceLength]
            }
        }
        finally {
            reader.close()
        }
    }
    
    void movingWindow(int windowSize, String chr, Closure c, Closure filterFn=null) {
        int chrSize = getContigs()[chr]
        movingWindow(windowSize, chr, 0, chrSize, c, filterFn)
    }
    
    
    /**
     * Call the given closure for each base position with a moving window of reads over that position
     * <p>
     * <b>NOTE:</b> Requires a BAM file sorted by position. Will not work with unsorted bam file.
     * 
     * @param windowSize    the size of window to use in bp
     * @param start         start position
     * @param end           end position
     * @param c             callback function to invoke
     * @param filterFn      optional filter function to apply that filters out reads
     */
    @CompileStatic
    void movingWindow(int windowSize, String chr, int start, int end, Closure c, Closure filterFn=null) {
        
       assert windowSize > 0
       
       ReadWindow readWindow = new ReadWindow()
        
       int halfWindowSize = (int)Math.floor((float)windowSize / 2.0f)
       
       int trailingPosition = start - halfWindowSize
       int leadingPosition = start + halfWindowSize
       int currentPosition = trailingPosition
       
       
       TreeMap<Integer,List<SAMRecord>> window = readWindow.window
       
       SAMFileReader reader = newReader()
       SAMRecordIterator iter = reader.query(chr, Math.max(0,start-windowSize), end+windowSize, false)
       try {
           try {
               while(iter.hasNext()) {
                   
                   SAMRecord r = iter.next()
                   
                   if((filterFn != null) && (filterFn(r) == false)) {
                       continue
                   }
                   
                   int pos = r.alignmentStart
                   
                   boolean verbose = false
                   if(pos == 18057685) {
                       println "Leading edge at debug position"
                   }
                   
                   if(currentPosition == 18057685) {
                       println "Window center at debug position"
                   }
     
                   
                   // Check if we already have reads at this position
                   List<SAMRecord> readsAtPos = window[pos]
                       
                   if(readsAtPos == null) { // encounter a new position, move window forward
                           
                       while(currentPosition<pos-halfWindowSize) {
                           
                           // Call for current position
                           readWindow.pos = currentPosition
                           if(currentPosition>=start && currentPosition<end)
                               c(readWindow)
                               
                           // Remove all that now fall outside lower window
                           int trailingEdge = currentPosition - halfWindowSize
                           while(!window.isEmpty() && window.firstKey()<trailingEdge) {
                               window.pollFirstEntry()
                           } 
                               
                           // Move forward 1
                           ++currentPosition
                       }
                       readsAtPos = window[pos] = []
                   }
                   readsAtPos << r
               }
               
               // We have ended where the reads ended, not necessarily where the requested
               // window ended. Keep invoking callback until we do it for every position 
               // requested
               while(currentPosition<end) {
                   
                   int trailingEdge = currentPosition - halfWindowSize
                   while(!window.isEmpty() && window.firstKey()<trailingEdge) {
                       window.pollFirstEntry()
                   }
                   
                   if(currentPosition == 18057685) {
                       println "Window center at debug position: reads are " + window[currentPosition]
                   }
      
                   readWindow.pos = currentPosition
                   if(currentPosition>=start && currentPosition<end)
                       c(readWindow)
                   ++currentPosition
               }
           }
           finally {
               iter.close()
           }
       }
       finally {
           reader.close()
       }
    }

    /**
     * Close the underlying SAMFileReader
     */
    void close() {
        if(this.samFileReader != null)
            this.samFileReader.close()
    }

    /**
     * Read a BAM or SAM file from standard input and call the given closure for
     * each read contained therein.
     * @param c
     */
    @CompileStatic
    static void eachRead(Closure c) {
        SAMFileReader reader = new SAMFileReader(System.in)
        Iterator<SAMRecord> iter = reader.iterator()
        while(iter.hasNext()) {
            c(iter.next())
        }
    }

    /**
     * Count the total number of reads in the SAM file
     * 
     * @return total number of reads (not pairs, each pair is counted as 2 reads)
     */
    @CompileStatic
    int size() {
        int count = 0
        eachRecord { SAMRecord r ->
            ++count
        }
        return count
    }
    
    /**
     * Probe the given BAM file to make a guess about what genome build it is 
     * generated from.
     * 
     * @return
     */
    String sniffGenomeBuild() {
        
        Map<String,Integer> contigs = getContigs()
        
        Map hgMap = [
            247249719 : "hg18",
            249250621 : "hg19",
            248956422 : "hg38"
        ]
        
        Map grcMap = [
            249250621 : "GRCh37",
            248956422 : "GRCh38"
        ]
        
        contigs.any { it.key.startsWith('chr') } ? 
            hgMap[contigs['chr1']]
        :
            grcMap[contigs['1']]
    }
    
    /*
    public static long getOpenFileDescriptorCount() {
        OperatingSystemMXBean osStats = ManagementFactory.getOperatingSystemMXBean();
        if(osStats instanceof UnixOperatingSystemMXBean) {
            return ((UnixOperatingSystemMXBean)osStats).getOpenFileDescriptorCount();
        }
        return 0;
    }
    */
    
    /**
     * Create an index for the given BAM file
     * 
     * @param bamFile
     */
    @CompileStatic
    static void index(File bamFile) {
            
        SAMFileReader reader = new SAMFileReader(bamFile)
        File outputFile = new File(bamFile.path + ".bai")
        
        if(outputFile.exists())
            outputFile.delete()
        
        BAMIndexer indexer = new BAMIndexer(outputFile, reader.getFileHeader());
            
        reader.enableFileSource(true);
        int totalRecords = 0;
            
        // create and write the content
        for (SAMRecord rec : reader) {
            indexer.processAlignment(rec);
        }
        indexer.finish();
    }
}
