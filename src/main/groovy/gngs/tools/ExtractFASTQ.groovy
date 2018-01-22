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
    
    Writer out
    
    public ExtractFASTQ(SAM bam, Regions regions) {
        super();
        this.bam = bam;
        this.regions = regions;
    }
    
    @CompileStatic
    void run(Writer out) {
        StringBuilder b = new StringBuilder(500)
        OrderedPairReader opr = new OrderedPairReader(this.bam)
        opr.progress.log = log
        opr.eachPair(regions) { SAMRecord r1, SAMRecord r2 ->
            b.setLength(0)
            
            // R1
            b.append("@${r1.readName} 1:N:0:1\n")
            b.append(r1.readString)
            b.append('\n+\n')
            b.append(r1.baseQualityString)
            
            // R2
            b.append("\n@${r1.readName} 2:N:0:1\n")
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
            gene 'Extract reads for the given gene symbol', args:Cli.UNLIMITED, required:false
            bam 'BAM file to extract reads from', args:1, required:true
        }
        
        OptionAccessor opts = cli.parse(args)
        if(!opts) {
            System.exit(1)
        }
        
        log.info "ExtractFASTQ starting ... "
        log.info "BAM=$opts.bam"
        
        Regions bed = null
        if(opts.L) {
            bed = new BED(opts.L).load()
            log.info "BED file $opts.L includes ${bed.size()} bp to scan"
        }
        else
        if(opts.genes) {
            String genomeVersion = new SAM(opts.bam).sniffGenomeBuild()
            RefGenes refGenes = RefGenes.download(genomeVersion)
            bed = new Regions()
            for(String gene in opts.genes) {
                Region geneRegion = refGenes.getGeneRegion(gene)
                if(geneRegion == null) {
                    System.err.println "Gene $opts.gene could not be resolved from RefSeq!"
                    System.exit(1)
                }
                log.info "Gene symbol $opts.gene translated to $geneRegion"
                bed.addRegion(geneRegion)
            }
        }
        
        System.out.withWriter { Writer w ->
            new ExtractFASTQ(new SAM(opts.bam), bed).run(w)
        }
    }
}
