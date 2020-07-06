package gngs.tools

import gngs.*
import graxxia.IntegerStats
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import groovyx.gpars.actor.Actor
import groovyx.gpars.actor.Actors
import htsjdk.samtools.BAMIndex
import htsjdk.samtools.BAMIndexMetaData
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReader
import htsjdk.samtools.SamReader.Indexing

@CompileStatic
class ReadSpans {
    int [][] intervals
    int count
    
    String toString() { "$count read spans from ${intervals[0][0]} to ${intervals[count-1][1]}"}
}

/**
 * SingleCov is a counterpoint to MultiCov, offering a subset of its functionality
 * for single samples at higher performance.
 * 
 * SingleCov produces per-base coverage and overall coverage statistics for a set 
 * of target regions on a BAM file. It writes a sample_summary file in the same
 * format that GATKDepthOfCoverage does, and per-base coverage in plain TSV,
 * GZIP or BGZIPPED format (auto-detected from the file name).
 * 
 * The main reason SingleCov is faster is because it trades memory for speed; it 
 * reads a whole chromosome of read spans at a time into an int-array and 
 * computes coverage depth from their overlaps in a block. This avoids all 
 * latency in reading, etc. This makes it 2-3 times faster than MultiCov for
 * single samples. However it takes several GB in memory to store the full
 * matrix of read spans and coverage values for each chromosome in memory,
 * making it a non-viable approach for multiple samples.
 * 
 * @author Simon Sadedin
 */
@Log
class SingleCov extends ToolBase {

    Regions scanRegions
    
    List<String> contigs
    
    Map<String,Integer> contigSizes
    
    IntegerStats covStats = new IntegerStats(1000)
    
    SAM bam 
    
    Writer out
    
    int minimumMappingQuality = 1
    
    @Override
    public void run() {
        
        if(!opts.arguments()) 
            throw new IllegalArgumentException("Please provide a BAM file to analyse")
            
        Regions targetRegions = new BED(opts.L).load()
        
        this.scanRegions = targetRegions.reduce().toSorted(new NumericRegionComparator())    
        
        this.contigs = scanRegions.index*.key
        
        this.bam = new SAM(opts.arguments()[0])
        
        this.contigSizes = bam.contigs
        
        if(opts.minMQ)
            this.minimumMappingQuality = opts.minMQ.toInteger()
            
        log.info "Analysing ${Utils.humanBp(scanRegions.size())} from ${opts.arguments()[0]}"
        log.info "Mapping quality threshold = $minimumMappingQuality"
        
        out = new AsciiWriter(Utils.outputStream(opts.o,10*1024*1024), 5*1024*1024)
        
        RegulatingActor writerActor = RegulatingActor.actor { contig, covs ->
            writeCoverage(contig, covs)
        }
        writerActor.start()
        
        RegulatingActor countActor = RegulatingActor.actor { contig, spans ->
            def covs = countCoverage(spans)
            writerActor.sendTo([contig,covs])
        }
        countActor.start()
        
        try {
            for(String contig in contigs) {
                log.info "Process $contig"
                ReadSpans spans = readContig(contig)
                countActor.sendTo([contig, spans])
            }

            countActor.sendStop()
            countActor.join()

            writerActor.sendStop()
            writerActor.join()
        }
        finally {
            out.close()
        }
        
        log.info "Coverage statistics: $covStats"
        
        if(opts.samplesummary) {
            writeSampleSummary(opts.samplesummary, covStats)
        }
    }
    
   
    @CompileStatic
    void writeCoverage(String contig, short [] covs) {
        final Regions contigRegions = scanRegions.getContigRegions(contig)
        log.info "Write ${Utils.human(contigRegions.size())} coverage values for contig: $contig"
        for(Region r in contigRegions) {
            final int start = r.from
            final int end = r.to
            final int zero = 0
            final int covLen = covs.size()
            for(int pos = start; pos < end; ++pos) {

                int cov = pos< covLen ? Math.max(Math.min(1000,covs[pos]),zero) : 0

                covStats.addIntValue(cov)
                out.write(contig)
                out.write('\t')
                out.write(String.valueOf(pos))
                out.write('\t')
                out.write(String.valueOf(cov))
                out.write('\n')
            }
        }
    }
    
    @CompileStatic
    short [] countCoverage(final ReadSpans spans) {

        int covSize = spans.intervals[spans.count-1][1]
        
        log.info "Allocate ${Utils.human(covSize*2)} bytes for coverage values"

        short [] covs =  new short[covSize]
        
        final ot = new gngs.OverlapTracker()
    
        int pos = 0
        
        final int[][] intervals = spans.intervals
        
        final int len = spans.count
    
        for(int i=0; i<len; ++i) {
    
            final int start = intervals[i][0]
            while(pos < start) {
                ot.removeNonOverlaps(pos)
                covs[pos] = (short)ot.size()
                ++pos
            }
    
            ot.add(intervals[i])
        }
        log.info "Computed ${covs.size()} coverage values"
        return covs
    }
    
    @CompileStatic
    ReadSpans readContig(final String contig) {
        
        int i = 0
        
        int recordCount = (int)bam.withReader { SamReader r ->
            
            def seqIndex = r.fileHeader.sequenceDictionary.sequences.find { it.sequenceName == contig }.sequenceIndex

            BAMIndexMetaData meta = ((Indexing)r).index.getMetaData(seqIndex)
            
            meta.getAlignedRecordCount()
        }

        def intervals = new int[recordCount][] as int[][]

        Region region = new Region("$contig:0-${contigSizes[contig]}")
        bam.withIterator(region) {
            while(it.hasNext()) {
                SAMRecord r = it.next()
                if(r.getReadUnmappedFlag())
                    continue
                if(r.secondaryOrSupplementary)
                    continue
                if(r.mappingQuality<minimumMappingQuality)
                    continue
                    
                final int mateStart = r.getMateAlignmentStart();
                int alignmentEnd = r.alignmentEnd
                int alignmentStart = r.alignmentStart

                // Don't double count if reads overlap
                if(r.getFirstOfPairFlag() && mateStart == alignmentStart) { 
                    // Note this logic is special: for the special case where the reads exactly overlap,
                    // there is no reliable choice which read is "first", so we can't decide which one
                    // to clip based on start position. That is why for the special case that they exactly
                    // overlap, we arbitrarily choose to exclude the first-of-pair completely from the coverage
                    // count and keep the whole second read
                    continue
                }
                else
                if(mateStart > alignmentStart && mateStart <= alignmentEnd) 
                    alignmentEnd = mateStart

                intervals[i] = [alignmentStart, alignmentEnd] as int[]
                ++i
            }
        }
        log.info "Scanned ${Utils.human(i)} reads from $bam.samFile"
        return new ReadSpans(intervals:intervals, count:i)
    }
    
    void writeSampleSummary(String outputFile, IntegerStats stats) {

        List<String> headers = ['Median Coverage', 'Mean Coverage', 'perc_bases_above_1', 'perc_bases_above_10', 'perc_bases_above_20', 'perc_bases_above_50']

        new File(outputFile).withWriter { w ->
            w.println(headers.join("\t"))
            w.println([
                stats.median,
                Utils.humanNumberFormat.format(stats.mean),
                Utils.humanNumberFormat.format(100*stats.fractionAbove(1)),
                Utils.humanNumberFormat.format(100*stats.fractionAbove(10)),
                Utils.humanNumberFormat.format(100*stats.fractionAbove(20)),
                Utils.humanNumberFormat.format(100*stats.fractionAbove(50)),
            ].join('\t'))
        }
        log.info "Wrote sample summary for ${bam.samples[0]} to $outputFile"
    }
    
    
    static void main(String[] args) {
        cli('SingleCov [-o] -L <target regions> <bam file>', args) {
            o 'Output file to write to', args:1, required: true
            minMQ 'Minimum mapping quality (1)', args:1, required: false
            samplesummary 'File to write coverage statistics to in tab separated format', args:1, required: false
            'L' 'bam file to read from', args:1, required: true, longOpt: 'target'
        }
    }
    
    
}
