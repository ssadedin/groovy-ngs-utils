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
package gngs.tools

import gngs.*
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import htsjdk.samtools.SAMRecord

/**
 * Extract reads as pairs from a sorted BAM and writes them out to standard 
 * output in interleaved format.
 * <p>
 * Example:
 * <p>
 * <code>gngstool ExtractFASTQ -bam my.bam -L chr1:5000000-6000000 > interleaved.fastq</code>
 * 
 * An example of realigning a portion of a BAM file:
 * <code>
 * gngstool ExtractFASTQ -bam my.bam -L chr1:5000000-6000000 | \
 *          bwa mem -p ref.fa  - | \
 *          samtools view -bSu - | \
 *          samtools sort -o out.bam 
 * </code>
 * 
 * @author Simon Sadedin
 */
@Log
class ExtractFASTQ {
    
    SAM bam
    
    Regions regions
    
    boolean addPositionToNames = false
    
    int spoolSize = 32000i
    
    int concurrency = 10
    
    Writer out
    
    public ExtractFASTQ(SAM bam, Regions regions) {
        super();
        this.bam = bam;
        this.regions = regions;
    }
    
    @CompileStatic
    void run(Writer out) {
        
        log.info "Using spool size $spoolSize, concurrency $concurrency"
        
        StringBuilder b = new StringBuilder(500)
        OrderedPairReader opr = new OrderedPairReader(this.bam, spoolSize:spoolSize, lookupConcurrency:concurrency)
        opr.progress.log = log
        opr.includeUnmapped = true
        opr.includeChimeric = true
        opr.eachPair(regions) { SAMRecord r1, SAMRecord r2 ->
            b.setLength(0)
            
            String r1Name = r1.readName
            if(addPositionToNames) {
                r1Name = r1Name + ":$r1.referenceName:$r1.alignmentStart:$r2.referenceName:$r2.alignmentStart"
            }
            
            // R1
            b.append("@${r1Name} 1:N:0:1\n")
            b.append(r1.readString)
            b.append('\n+\n')
            b.append(r1.baseQualityString)
            
            // R2
            b.append("\n@${r1Name} 2:N:0:1\n")
            b.append(gngs.FASTA.reverseComplement(r2.readString))
            b.append('\n+\n')
            b.append(r2.baseQualityString.reverse())
            
            out.println(b.toString())
        }
    }
    
    static void main(String [] args) {
        
        Utils.configureSimpleLogging()
        
        Cli cli = new Cli(usage: 'groovy gngs.tools.ExtractFASTQ -bam <bam>')
        cli.with {
            'L' 'Region to extract', args:1, required:false
            pad 'Amount to pad regions by (0)', args:1, required: false
            gene 'Extract reads for the given gene symbol', args:Cli.UNLIMITED, required:false
            buffer 'Size of buffer to use for pairing reads (32000)', args:1, required:false
            namepos 'Add original position to the read names', required:false
            bam 'BAM file to extract reads from', args:1, required:true
            o 'Output file (default = stdout)', longOpt: 'output', args:1, required:false
        }
        
        OptionAccessor opts = cli.parse(args)
        if(!opts) {
            System.exit(1)
        }
        
        log.info "ExtractFASTQ starting ... "
        log.info "BAM=$opts.bam"
        
        Regions bed = null
        if(opts.L) {
            bed = new BED(opts.L).load().reduce()
            log.info "BED file $opts.L includes ${bed.size()}bp (${Utils.humanBp(bed.size())}) to scan"
        }
        else
        if(opts.genes) {
            String genomeVersion = new SAM(opts.bam).sniffGenomeBuild()
            bed = getGeneRegions(genomeVersion, opts.genes)
        }
        
        if(bed != null && opts.pad) {
            log.info "Padding regions by ${opts.pad}bp"
            bed = bed.widen(opts.pad.toInteger())
            bed = bed.reduce()
            log.info "After padding, regions span ${bed.size()}bp (${Utils.humanBp(bed.size())})"
        }
        
        Writer outputWriter = opts.o ? new File(opts.o).newWriter() : System.out.newWriter()
        outputWriter.withWriter { Writer w ->
            ExtractFASTQ efq = new ExtractFASTQ(new SAM(opts.bam), bed)
            if(opts.namepos)
                efq.addPositionToNames = true
            if(opts.buffer) 
                efq.spoolSize = opts.buffer.toInteger()
            efq.run(w)
        }
    }
    
    /**
     * Attempt to interpret the given list of genes as either a list of gene symbols
     * or a text file specifying the genes to read.
     * 
     * @param genomeVersion
     * @param genes
     * @return
     */
    static Regions getGeneRegions(String genomeVersion, List<String> genes) {
        RefGenes refGenes = RefGenes.download(genomeVersion)
        Regions geneRegions = new Regions()
        List<String> geneList = genes
        if(new File(geneList).exists()) {
            log.info "$geneList exists as file: interpreting as file to read gene symbols from"
            geneList = new File(geneList).readLines()*.trim()
            log.info "Read ${geneList.size()} gene symbols"
        }
            
        for(String gene in geneList) {
            Region geneRegion = refGenes.getGeneRegion(gene)
            if(geneRegion == null) {
                log.warning "Gene $gene could not be resolved from RefSeq!" 
            }
            log.info "Gene symbol $gene translated to $geneRegion"
            geneRegions.addRegion(geneRegion)
        }
        return geneRegions
    }
}
