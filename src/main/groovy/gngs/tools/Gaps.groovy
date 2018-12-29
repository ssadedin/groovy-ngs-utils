/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2018 Simon Sadedin, ssadedin<at>gmail.com
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
package gngs.tools

import java.text.NumberFormat

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics

import gngs.*
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import groovyx.gpars.GParsPool
import groovyx.gpars.GParsPoolUtil

@Log
class GapAnnotator extends RegulatingActor<CoverageBlock> {
    
    RefGenes refgenes

    public GapAnnotator(RefGenes refgenes) {
        super(1000, 10000);
        this.refgenes = refgenes
    }
    
    final static List<String> ANNOTATION_COLUMNS = ['transcript', 'gene', 'exon', 'coding_intersect']

    @Override
    public void process(CoverageBlock block) {
        
        // Check for overlapping transcripts
        List<GRange> transcripts = refgenes.refData.getOverlaps(block)
        block.annotations = []
        if(transcripts) {
            for(Region info in transcripts*.extra) {
                block.id = info.gene
                
                int cdsStart = info.cds_start.toInteger()
                int cdsEnd = info.cds_end.toInteger()
                Region cdsRegion = new Region(block.chr, cdsStart, cdsEnd)
                Regions exons = refgenes.getTranscriptExons(info.tx)
                Regions intersectedExons = exons.intersect(new Regions([block]))
//                log.info "intersected Exon numbers: " + intersectedExons*.exon
                for(Region exon in exons) {
                    if(!exon.overlaps(block))
                        continue
                        
                    int cdsIntersect = Math.max(exon.intersect(cdsRegion).intersect(block).size()-1,0)
                    
                    block.annotations << [ 
                        transcript: info.tx, 
                        gene: info.gene, 
                        exon: info.strand == "+" ? exon.exon : (exons.numberOfRanges - exon.exon),
                        coding_intersect: cdsIntersect,
                    ]
                }
                
            }
//            log.info "Annotated gap: " + block
        }
        
    }
}

/**
 * A tool that computes gaps in coverage based on BEDTools or MultiCov format output.
 * <p>
 * Can compute simple gap calculations, but also diff gaps between two samples
 * and also between sets of samples.
 * <p>
 * Simple usage:
 * <pre>
 * gngstool Gaps -h coverage.txt
 * </pre>
 * Diff usage:
 * <pre>
 * gngstool Gaps -h -diff oldcoverage1.txt,oldcoverag2 coverage1.txt coverage2.txt
 * </pre>
 * 
 * @author Simon Sadedin
 */
@Log
class Gaps {
    
    OptionAccessor opts
    
    Regions targetRegions
    
    int concurrency = 2
    
    NumberFormat format = NumberFormat.numberInstance
    
    Gaps(OptionAccessor opts) {
        this.opts = opts
        this.targetRegions = null
        if(opts.L) {
            this.targetRegions = new BED(opts.L,withExtra:true).load(withExtra:true).collect {
                new Region(it.chr, it.from+1, it.to, id: it.range.extra)
            } as Regions
        }
        if(opts.n)
            concurrency = opts.n.toInteger()
            
        format.maximumFractionDigits = 1
        format.minimumFractionDigits = 1
    }
    
    void run() {
        List<Regions> allGaps = computeAllGaps()
        
        Regions unionGaps = computeGapUnion(allGaps)
        
        if(opts.diff) {
            List<Region> oldGaps = computeDiffGaps()
            
            Regions unionDiffGaps = computeGapUnion(oldGaps)
            writeDiffs(unionDiffGaps, unionGaps)
        }
        else {
            log.info "Writing output ..."
            if(opts.h)
                writeTable(unionGaps)
            else {
                writeGaps(unionGaps)
                log.info "Done."
            }
        }
    }
    
    @CompileStatic
    Regions computeGapUnion(List<Regions> allGaps) {
        
        Regions unionGaps
        if(allGaps.size() == 1)
            unionGaps = allGaps[0]
        else {
            unionGaps = new Regions()
            for(Regions gaps in allGaps) {
                for(Region r in gaps) {
                    unionGaps.addRegion(r)
                }
            }
        }
        
        Regions reduced = unionGaps.reduce()
        
        Regions filtered = reduced
        if(this.targetRegions != null) {
            filtered = new Regions()
            for(Region r in reduced) {
                for(Range ix in targetRegions.intersect(r)) {
                    if(ix.from == ix.to) {
                        continue
                    }
                    
                    Region ixRegion = (Region)((GRange)ix).extra
                    filtered.addRegion(new Region(
                        r.chr, 
                        new GRange((int)ix.from, (int)ix.to, 
                            new CoverageBlock(chr: r.chr, start: (int)ix.from, end: (int)ix.to, 
                                stats: ((DescriptiveStatistics)((CoverageBlock)((GRange)r.range).extra).stats), 
                                id: (String)ixRegion.getProperty('id'))
                        )
                    ))
                }
            }
        }
        
        return filtered
    }
    
    List<Regions> computeAllGaps() {
        return GParsPool.withPool(concurrency) { 
            opts.arguments().collectParallel { String coverageFile ->
                log.info "Calculating coverage gaps for $coverageFile ..."
                CoverageGaps gaps = calculateGaps(coverageFile) 
                return gaps as Regions 
            }
        }
    }
    
    List<Regions> computeDiffGaps() {
        return GParsPool.withPool(concurrency) { 
            return opts.diff.tokenize(',')*.trim().collectParallel { String coverageFile ->
                log.info "Calculating coverage gaps for $coverageFile ..."
                CoverageGaps gaps = calculateGaps(coverageFile) 
                return gaps as Regions
            }
        }
    } 
    
    void logGapInfo(String file, CoverageGaps gaps) {
        log.info "Median coverage for $file = $gaps.coveragePercentiles.median"
        log.info "Key percentiles: " + [20,30,50,100].collect { 
            format.format(100*gaps.coveragePercentiles.fractionAbove(it)) + ">" +it+"x" 
        }.join(', ')
    }
    
    public static void main(String [] args) {
        
        Utils.configureSimpleLogging()
        
        Cli cli = new Cli(usage: "Gaps <BEDTools coverage file>")
        cli.with {
            t 'Coverage threshold', longOpt: 'threshold', args:1, required:false
            diff 'Show only gaps not occurring in file(s)', args:1, required:false
            'L' 'Only output results for <bed file> regions', longOpt: 'regions', args:1, required:false
            n 'Concurrency to use', args:1, required:false
            m 'Input file is multicov format, not BEDTools'
            r 'Annotate with refgene for genome build <arg>', args:1
            a 'Write annotated report to separate file <arg>', args:1
            h 'Output a human readable table'
        }
        
        OptionAccessor opts = cli.parse(args)
        if(!opts.arguments()) {
            cli.usage()
            System.err.println "Please provide a BEDTools coverage file"
            System.exit(1)
        }
        
        Gaps gaps = new Gaps(opts)
        gaps.run()
    }
    
    CoverageGaps calculateGaps(String coverageFile) {
        CoverageGaps gaps = new CoverageGaps(coverageFile)
        if(opts.r) {
            RefGenes refgenes = new File(opts.r).exists() ? new RefGenes(opts.r) : RefGenes.download(opts.r)
            gaps.gapProcessor = new GapAnnotator(refgenes)
            gaps.gapProcessor.start()
        }
        
        if(opts.t)
            gaps.threshold = opts.t.toInteger()
            
        if(opts.m)
            gaps.calculateMultiCov()
        else
            gaps.calculate()
            
        logGapInfo(coverageFile,gaps)
        return gaps
    }
    
    void writeDiffs(Regions oldRegions, Regions newRegions) {
        
        Regions introducedGaps = newRegions.subtract(oldRegions)
        introducedGaps = introducedGaps.enhance().grep { it.to > it.from } as Regions
        introducedGaps.each { it.type = 'introduced' }
        
        Regions eliminatedGaps = oldRegions.subtract(newRegions).grep { it.to > it.from } as Regions
        eliminatedGaps = eliminatedGaps.enhance()
        eliminatedGaps.each { it.type = 'eliminated' }
        
        Regions diffs = introducedGaps + eliminatedGaps
        log.info "Detected ${diffs.numberOfRanges} different gaps ($introducedGaps.numberOfRanges introduced, $eliminatedGaps.numberOfRanges eliminated)"
        log.info "Size of differences: ${diffs.size()}, Size of introduced: ${introducedGaps.size()}, Size of eliminated: ${eliminatedGaps.size()}"
        
        if(opts.h) {
            Utils.table(
               ["Chr","Start","End","Length","Difference"],
               diffs.collect {  Region block ->
                   [block.chr, block.from, block.to, block.size(), block.type] 
               }
        )}
        else {
            for(Region diff in diffs) {
                println([diff.chr, diff.from, diff.to, diff.size(), diff.type].join('\t'))
            }
        }
    }
    
    void writeTable(Regions gaps) {
        Utils.table(
           ["Chr","Start","End","ID","Length","Min","Mean","Max"],
           gaps.collect {  Region region ->
               CoverageBlock block = region.extra
               [block.chr, block.start, block.end, block.id, block.end - block.start, (int)block.stats.min, block.stats.mean, (int)block.stats.max] 
           }
        )
    }
    
    void writeGaps(Regions gaps) {
        
        Writer annotationWriter
        if(opts.a) {
            annotationWriter = Utils.writer(opts.a)
        }
        
        
        try {
            for(Region blockRegion in gaps) {
                writeGap(blockRegion, annotationWriter)
            }
        }
        finally {
            if(annotationWriter != null)
                annotationWriter.close()
        }
    }
    
    /**
     * Write the given coverage gap to the output.
     * <p>
     * If the gap is annotated, and a separate annotationWriter is provided, the annotated
     * form will be written to the annotationWriter. Otherwise, annotations (if they exist)
     * are written to standard output along with the gaps.
     * 
     * @param blockRegion
     * @param annotationWriter
     */
//    @CompileStatic
    void writeGap(Region blockRegion, Writer annotationWriter) {
        CoverageBlock block = blockRegion.extra
        List fields = [block.chr, block.start, block.end, block.id, block.end - block.start, (int)block.stats.min, block.stats.mean, (int)block.stats.max]
        if(opts.r) {
            if(block.annotations) {
                for(Map annotation in block.annotations) {
                    List rowFields = fields + GapAnnotator.ANNOTATION_COLUMNS.collect { annotation[it] }
                    if(annotationWriter != null) {
                        annotationWriter.println(rowFields.join('\t'))
                        println(fields.join('\t'))
                    }
                    else
                        println(rowFields.join('\t'))
                }
            }
            else { // Add empty values for the annotations
                String annotatedFields = (fields + (['']*GapAnnotator.ANNOTATION_COLUMNS.size())).join('\t')
                if(annotationWriter != null) {
                    annotationWriter.println(annotatedFields)
                    println(fields.join('\t'))
                }
                else {
                    println(annotatedFields)
                }
            }
        }
        else {
            println(fields.join('\t'))
        }
    }

}
