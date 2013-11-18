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

import com.sun.management.UnixOperatingSystemMXBean;

import groovy.transform.CompileStatic;
import net.sf.samtools.BAMRecord;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;


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
}

/**
 * Utility class for handling SAM / BAM files
 * 
 * @author simon.sadedin@mcri.edu.au
 *
 */
class SAM {
    
    SAMFileReader samFileReader;
    
    File samFile
    
    File indexFile
    
    int minMappingQuality = 1
    
    SAM(String fileName) {
        this(new File(fileName))
    }
    
    SAM(File file) {
        this.samFile = file
        if(!this.samFile.exists())
            throw new FileNotFoundException("BAM file could not be opened at ${samFile.absolutePath}")
        this.indexFile = new File(samFile.absolutePath + ".bai")
        if(!indexFile.exists()) {
            indexFile = new File(samFile.absolutePath.replaceAll(".bam\$",".bai"))
        }
        if(!indexFile.exists())
            throw new FileNotFoundException("Please ensure your BAM / SAM file is indexed. File $indexFile could not be found.")
        
        SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT)
        this.samFileReader = new SAMFileReader(samFile, indexFile, false) 
    }
    
    /**
     * Iterate over each record in the same file in the order they are in the file
     * @param c
     */
    void eachRecord(Closure c) {
        use(SAMRecordCategory) {
          Iterator i = samFileReader.iterator()
          while(i.hasNext()) {
              c(i.next())
          }
        }
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
    DescriptiveStatistics coverageStatistics(String chr, int pos, int end) {
        DescriptiveStatistics stats = new DescriptiveStatistics()
        int total = 0
        this.pileup(chr, pos, end) { PileupIterator.Pileup p ->
            stats.addValue(p.alignments.size())
        }
        stats
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
        SAMRecordIterator i = r.query(chr, pos, end, false);
        try {
            int count = 0;
            while(i.hasNext()) {
                SAMRecord rec = (SAMRecord)i.next();
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
    
    void close() {
        if(this.samFileReader != null)
            this.samFileReader.close()
    }
    
    @CompileStatic
    static void eachRead(Closure c) {
        SAMFileReader reader = new SAMFileReader(System.in)
        Iterator<SAMRecord> iter = reader.iterator()
        while(iter.hasNext()) {
            c(iter.next())
        }
    }
    
    public static long getOpenFileDescriptorCount() {
        OperatingSystemMXBean osStats = ManagementFactory.getOperatingSystemMXBean();
        if(osStats instanceof UnixOperatingSystemMXBean) {
           return ((UnixOperatingSystemMXBean)osStats).getOpenFileDescriptorCount();
        }
        return 0;
    }
}
