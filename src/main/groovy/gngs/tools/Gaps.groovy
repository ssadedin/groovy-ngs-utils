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
import groovy.transform.Memoized
import groovy.util.logging.Log
import groovyx.gpars.GParsPool
import groovyx.gpars.GParsPoolUtil
import graxxia.TSV

/**
 * Annotates regions with gene, transcript, CDS overlap and distance information, as well
 * as custom annotation to gene associations ("panels") and sub-panels 
 * (groups of genes within a panel)
 * 
 * @author Simon Sadedin
 */
@Log
class GapAnnotator extends RegulatingActor<CoverageBlock> {
    
    RefGenes refgenes

    /**
     * Map keyed on gene to a nested lookup table that is keyed on panel name and returns the list of 
     * subpanels from the major panel that include the gene.
     */
    Map<String, Map<String, List>> panelGeneMap = new HashMap<String, Map<String, List>>()

    /**
     * List of all known sub panels ("panel classes")
     */
    ArrayList<String> panelClasses = []


    public GapAnnotator(RefGenes refgenes) {
        super(1000, 10000);
        this.refgenes = refgenes
    }

    public GapAnnotator(RefGenes refgenes, ArrayList<String> panelFiles) {
        this(refgenes)

        indexPanels(panelFiles)
    }
    
    void indexPanels(List<String> panelFiles) {
        
        List panelClassNames = panelFiles.collect {

            def fileName = new File(it).name.trim()

            // File name will be the upperclass name for the panel
            fileName.take(fileName.lastIndexOf('.'))
        }
        
        indexPanels(panelClassNames, panelFiles.collect { new TSV(it).toListMap() })
    }
    
    void indexPanels(List<String> panelClassNames, List<List<Map>> tsvs) {

        // Process each panel File 
        [tsvs, panelClassNames].transpose().each { lm, panelClassName ->

            this.panelClasses.add(panelClassName)

            lm.each { Map row ->
                row.each { rawSubClassName, gene -> 
                    
                    def subClassName = rawSubClassName.trim()
                    
                    if(!gene) 
                        return

                    if (!this.panelGeneMap.containsKey(gene)) {
                        this.panelGeneMap[gene] = [:]
                    }
                    
                    List geneSubClasses = this.panelGeneMap[gene].computeIfAbsent(panelClassName) { [] }

                    if(!(subClassName in geneSubClasses))
                        geneSubClasses << subClassName
                }
            }
        }
        
    }
    
    /**
     * Names of properties set as annotations on gap regions
     */
    final static List<String> ANNOTATION_COLUMNS = ['transcript', 'strand', 'gene', 'exon', 'coding_intersect','cds_distance','transcript_type']
    
    /**
     * Prettier human readable version of above columns - should always match 1-1
     */
    final static List<String> ANNOTATION_OUTPUT_COLUMNS = ['Tx Name', 'Strand', 'Symbol', 'Exon Number', 'CDS Overlap', 'CDS Distance', 'Transcript Type']

    @CompileStatic
    @Override
    public void process(CoverageBlock block) {
        
        final Region blockRegion = new Region(block)
        
        // Check for overlapping transcripts
        assert block.annotations == null || block.annotations.isEmpty()
        block.annotations = this.annotateGapRegion(blockRegion)
        block.id = block.annotations[0]?.gene?:'N/A'
    }
    
    @CompileStatic
    List<Map> annotateGapRegion(final Region blockRegion) {

        List<GRange> transcripts = (List<GRange>)refgenes.refData.getOverlaps(blockRegion)

        if(transcripts == null || transcripts.isEmpty()) {
            return []
        }

        List<Map> annotations = []
            
        List<Region> transcriptRegions = new ArrayList(transcripts.size())
        for(GRange gr in transcripts) {
            transcriptRegions.add((Region)gr.extra)
        }
        
        for(final Region transcriptRegion in transcriptRegions) {
            
            Map info = transcriptRegion.properties
           
            int cdsStart = info['cds_start'].toString().toInteger()
            int cdsEnd = info.cds_end.toString().toInteger()

            Region cdsRegion = new Region(blockRegion.chr, cdsStart, cdsEnd)
            Regions cdsRegions = new Regions()
            cdsRegions.addRegion(cdsRegion)

            Regions allTranscriptExons = getExonsForTranscript((String)info.tx)
            Regions intersectedExons = allTranscriptExons.intersect(new Regions([(IRegion)blockRegion]))
            Regions codingExons = allTranscriptExons.intersectRegions(cdsRegions)

            List transcriptAnnotations = []

            Map panelAnnotations = getPanelAnnotations((String)info.gene)

            for(Region exon in allTranscriptExons) {
                if(!blockRegion.overlaps(exon))
                    continue

                def annotation = [ 
                    id: info.gene,
                    transcript: info.tx, 
                    strand: info.strand,
                    gene: info.gene, 
                    exon: info.strand == "+" ? exon['exon'] : (allTranscriptExons.numberOfRanges - (Integer)exon['exon']),
                    cds_distance: codingExons.distanceTo(blockRegion),
                    transcript_type: (cdsStart == cdsEnd) ? 'non-coding' : 'coding',
                    coding_intersect: exon.overlaps(blockRegion) ? Math.max((int)exon.intersect(cdsRegion).intersect(blockRegion).size()-1,0) : 0,
                ] + panelAnnotations

                transcriptAnnotations << annotation
            }
            
            // If no exon from the transcript actually overlapped, then we write out a single line annotation line for the transcript
            if(transcriptAnnotations.isEmpty()) {
                def transcriptAnnotation = [
                    id: info.gene,
                    transcript: info.tx, 
                    strand: info.strand,
                    gene: info.gene, 
                    exon: 'N/A',
                    cds_distance: codingExons.distanceTo(blockRegion),
                    transcript_type: (cdsStart == cdsEnd) ? 'non-coding' : 'coding',
                    coding_intersect: 0
                ] + panelAnnotations

                transcriptAnnotations << transcriptAnnotation
            }
            

            
            
            annotations.addAll(transcriptAnnotations)
        }

        // If any had coding sequence overlap then return only the set that have that
        // Otherwise return the annotation with minimum cds distance
        if(annotations.isEmpty()) {
            return []
        }
        else
        if(annotations.any { it.cds_distance == 0 }) {
            return annotations.grep { Map ann -> ann.cds_distance == 0 }
        }
        else {
            return [annotations.min { Map ann -> ann.cds_distance }]
        }
    }
    
    @CompileStatic
    @Memoized(maxCacheSize=500)
    Map<String,Object> getPanelAnnotations(String gene) {
        return panelClasses.collectEntries { panel ->
            if (!panelGeneMap[gene]?.containsKey(panel)) {
                [panel, 'no']
            } else {
                List subClasses = (List)panelGeneMap[gene][panel]
                [panel, subClasses.join(',')]
            }
        }
    }
    
    @CompileStatic
    @Memoized(maxCacheSize=500)
    Regions getExonsForTranscript(String tx) {
       refgenes.getTranscriptExons(tx) 
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
    
    final static List DEFAULT_COLUMNS = ['Chr', 'Start', 'End', 'Gene', 'Width', 'Min Cov', 'Mean Cov', 'Max Cov']
    
    CliOptions opts
    
    Regions targetRegions
    
    int concurrency = 2
    
    NumberFormat format = NumberFormat.numberInstance
    
    Writer gapWriter 
    
    Gaps(CliOptions opts) {
        this.opts = opts
        this.targetRegions = null
        this.gapWriter = new PrintWriter(System.out)
        if(opts.L) {
            this.targetRegions = new BED(opts.L,withExtra:true).load(withExtra:true).collect {
                new Region(it.chr, it.from+1, it.to, id: it.range.extra)
            } as Regions
        }
        if(opts.n)
            concurrency = opts.n.toInteger()
            
        format.maximumFractionDigits = 1
        format.minimumFractionDigits = 0
        format.groupingUsed = false
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
                    CoverageBlock block = (CoverageBlock)((GRange)r.range).extra
                    filtered.addRegion(new Region(
                        r.chr, 
                        new GRange((int)ix.from, (int)ix.to, 
                            new CoverageBlock(chr: r.chr, start: (int)ix.from, end: (int)ix.to, 
                                stats: ((DescriptiveStatistics)block.stats), 
                                id: (String)ixRegion.getProperty('id'),
                                annotations: block.annotations)
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
            csv 'Write output comma separated instead of tab separated'
            h 'Output a human readable table'
            p 'Set .tsv files [PanelClassName].tsv with Subclass names as headers', args:Cli.UNLIMITED, required:false
        }

        
        OptionAccessor opts = cli.parse(args)
        if(!opts) 
            System.exit(1)

        if(!opts.arguments()) {
            cli.usage()
            System.err.println "Please provide a BEDTools coverage file"
            System.exit(1)
        }

        if(opts.ps) {
            def err = false
            for(panelFile in opts.ps) {
                def currPanelFile = new File(panelFile)
                if (!currPanelFile.exists()) {
                    System.err.println "${panelFile} doesn't exist. Please provide an exisitng path"
                    err = true
                } 
            }
            if (err) {
                System.exit(1)
            }
        }
        
        Gaps gaps = new Gaps(new CliOptions(opts:opts))
        gaps.run()
    }
    
    CoverageGaps calculateGaps(String coverageFile) {
        CoverageGaps gaps = new CoverageGaps(coverageFile)
        if(opts.r) {
            RefGenes refgenes = new File(opts.r).exists() ? new RefGenes(opts.r) : RefGenes.download(opts.r)
            gaps.gapProcessor = opts.ps ? new GapAnnotator(refgenes, opts.ps) : new GapAnnotator(refgenes)
            gaps.gapProcessor.start()
            log.info "Annotating using gene definitions from $opts.r"
        }

        if(opts.t)
            gaps.threshold = opts.t.toInteger()
            
        
        if(opts.m)
            gaps.calculateMultiCov()
        else
            gaps.calculate()
            
//        logGapInfo(coverageFile,gaps)
        
        if(opts.r) {
            gaps.gapProcessor.sendStop()
            gaps.gapProcessor.join()
        }
            
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

        def panelClasses = null
        
        List cols = [] + DEFAULT_COLUMNS // clone
        if(opts.csv) {
            if(opts.r && !opts.a)
                cols += GapAnnotator.ANNOTATION_OUTPUT_COLUMNS
            if(opts.p) {
                panelClasses = opts.ps.collect {
                    def fileName = new File(it).name
                    fileName.trim().take(fileName.lastIndexOf('.'))
                }
                cols += panelClasses
            }

            println(cols.join(','))
        }
        
        try {
            for(Region blockRegion in gaps) {
                writeGap(blockRegion, annotationWriter, panelClasses)
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
     * @param panelClasses
     */
    @CompileStatic
    void writeGap(Region blockRegion, Writer annotationWriter, List panelClasses = null) {
        CoverageBlock block = (CoverageBlock)blockRegion.extra
        writeGapBlock(block, annotationWriter, panelClasses)
    }

    public void writeGapBlock(CoverageBlock block, Writer annotationWriter, List panelClasses = null) {

        List fields = [block.chr, block.start, block.end, block.id, block.end - block.start, (int)block.stats.min, block.stats.mean, (int)block.stats.max]

        String sep = opts.csv ? ',' : '\t'

        if(opts.r) {
            //            log.info "Write gap (${block.hashCode()}): " + block
            if(block.annotations) {
                for(Map annotation in block.annotations) {
                    List rowFields = (fields + GapAnnotator.ANNOTATION_COLUMNS.collect { annotation[it] }).collect {
                        Utils.formatIfNumber(format, it)
                    }
                    // Add panel annotations
                    if (panelClasses != null && annotation.containsKey(GapAnnotator.PANEL_ANNOTATION_KEY)) {
                        panelClasses.each { c ->
                            rowFields << annotation[GapAnnotator.PANEL_ANNOTATION_KEY][c].trim()
                        }
                    }
                    if(annotationWriter != null) {
                        annotationWriter.println(rowFields.join(sep))
                        if(gapWriter)
                            gapWriter.println(fields.join(sep))
                    }
                    else
                        gapWriter.println(rowFields.join(sep))
                }
            }
            else { // Add empty values for the annotations
                String annotatedFields = (fields + (['']*GapAnnotator.ANNOTATION_COLUMNS.size())).join(sep)
                if(annotationWriter != null) {
                    annotationWriter.println(annotatedFields)
                    if(gapWriter)
                        gapWriter.println(fields.join(sep))
                }
                else {
                    if(gapWriter)
                        gapWriter.println(annotatedFields)
                }
            }
        }
        else {
            gapWriter.println(fields.join(sep))
        }
    }

}
