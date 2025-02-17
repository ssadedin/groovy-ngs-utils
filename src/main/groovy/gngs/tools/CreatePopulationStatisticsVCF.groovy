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

import gngs.*
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Log
import htsjdk.variant.variantcontext.Genotype
import htsjdk.variant.variantcontext.GenotypesContext
import htsjdk.variant.variantcontext.VariantContext

import static gngs.Sex.*

import java.text.DecimalFormat

/**
 * A utility class to contain some functionality for computing allele numbers / counts
 * 
 * @author Simon Sadedin
 */
@CompileStatic
class AlleleNumber {
   
    static boolean isX(String chr) {
        chr.endsWith('X')
    }
    
    static boolean isY(String chr) {
        chr.endsWith('Y')
    }
    
    static int getAlleleNumber(gngs.Sex sex, String chr) {
        
        if(sex == MALE) {
            if(isX(chr) || isY(chr))
                return 1
            else
                return 2
        }
        else
        if(sex == FEMALE) {
            if(isY(chr))
                return 0
            else
                return 2
        }
        else
            assert false : "Only MALE and FEMALE sexes are supported"
    }
}

/**
 * A tool to create a "population statistics VCF" that represents the frequency of
 * variants from a given set of samples by entries in the info field.
 * <p>
 * The following population metrics are calculated:
 * <li>AN   -   the total number of possible alleles that could have been observed at the site
 * <li>AC   -   the total number of alleles that were observed at the site
 * <li>GTC  -   the total number of individuals that were sampled to generate AN and AC
 * <p>
 * The current implementation infers the sexes of the samples so they do not need to be specified.
 * However this does rely on having sufficient coverage of the autosomes and X chromosome to determine
 * sex. If this is not true, the sex of the samples must be provided with `-sex` or an error will 
 * be reported,
 * <p> 
 * <b>Limitations:</b>
 * The current implementation has a number of limitations:
 * 
 * <li>Only VCF information is considered, not gVCF, which is to say, an allele is counted as 
 *     unobserved at a site even if there was no coverage and therefore was not ascertained at 
 *     that site.
 * 
 * These limitations should be corrected to generate truly accurate counts.
 * 
 * @author Simon Sadedin
 */
@Log
class CreatePopulationStatisticsVCF extends ToolBase {
    
    int numSamples
    
    /**
     * The inferred sexes as a list
     */
    private List<gngs.Sex> sexes
    
    /**
     * The inferred sexes of the samples in the VCFs
     */
    Map<String,gngs.Sex> sampleSexes = [:]
    
    PrintWriter out = null

    @Override
    public void run() {
        
        if(opts.arguments().isEmpty()) {
            parser.usage()
            System.exit(0)
        }
        
        if(out != null) {
            // Use it as is
        }
        else
        if(opts.o) {
            out = Utils.writer(opts.o)
        }
        else {
            out = new PrintWriter(System.out)
        }
        
        out.withWriter {

            printVCFHeader()
            
            out.flush()
            
            inferSexes()
            
            // Since we end up with a sex inferred for every sample ...
            numSamples = sampleSexes.size()
            
            VariantContext lastVariant = null
            ProgressCounter progress = new ProgressCounter(withRate: true, withTime:true, extra: {
                "chr: $lastVariant.contig pos: $lastVariant.start (${lastVariant.genotypes[0].sampleName})"
            })
            VCFSiteWalker locusWalker = new VCFSiteWalker(opts.arguments())
            locusWalker.walk { List<VariantContext> variants ->
                lastVariant = variants[0]
                processSite(variants)
                progress.count()
            }
            progress.end()
            out.close()
        }
    }
    
    @CompileStatic
    void processSite(List<VariantContext> variants) {
        // TODO: in here we should group the variants by their alternate allele and output
        // a different count for each allele (to deal with multi-allelic variants)
        
        Map<String, List<VariantContext>> alleles = 
            variants.groupBy { VariantContext ctx -> 
                ctx.getAlternateAllele(0).baseString
            }
            
        if(alleles.size()>1) {
            log.info "Site ${variants[0].contig}:${variants[0].start} is multiallelic"
        }
        
        for(List<VariantContext> allele : alleles.values()) {
            VariantContext v0 = allele[0]
            PopulationAlleleCounts pac = computeAlleleCount(allele)
            int an = computeAlleleNumber(v0.contig)
            int gtc = numSamples // todo: check this is right
            printVCFSite(v0, pac, an, gtc)
        }
    }
    
    @CompileStatic
    @Memoized
    int computeAlleleNumber(String contig) {
       (int)sexes.sum { gngs.Sex sex ->  AlleleNumber.getAlleleNumber(sex, contig) } 
    }
    
    DecimalFormat afFormat = new java.text.DecimalFormat("#.####")

    /**
     * Print out a line of the output VCF based on the given computed statistics
     * 
     * @param v0
     * @param ac
     * @param an
     */
    @CompileStatic
    protected void printVCFSite(VariantContext v0, PopulationAlleleCounts ac, int an, int gtc) {
        
        def af = an > 0 ? afFormat.format(ac.ac / an.toDouble()) : '0'
        
        out.println([
            v0.contig,
            v0.start,
            v0.getID(),
            v0.reference.baseString,
            v0.alternateAlleles[0],
            '.',
            '.',
            "AC=$ac.ac;AN=$an;GTC=$gtc;AF=$af;HOM=$ac.hom;HET=$ac.het;HEMI=0"
        ].join('\t'))
    }
    
    /**
     * Compute the allele count at a site, based on the given list of variants at that site.
     * 
     * @return  the number of alleles observed, based on the het / hom status of the variants
     *          and the sex of the respective samples
     */
    @CompileStatic
    PopulationAlleleCounts computeAlleleCount(List<VariantContext> variants) {
        Set<String> samples = new HashSet(numSamples*2)
        PopulationAlleleCounts ac = (PopulationAlleleCounts)variants.sum { VariantContext vc -> // sum across VCFs
            vc.genotypes.collect { Genotype gt ->         // sum across samples within this vcf
                
                // Make us robust to samples being provided twice
                if(samples.contains(gt.sampleName in samples))
                    return PopulationAlleleCounts.ZERO_COUNTS
                
                samples.add(gt.sampleName)
                
                gngs.Sex sex = sampleSexes[gt.sampleName]
                int sampleAN = AlleleNumber.getAlleleNumber(sex, vc.contig)
                GenotypesContext ctx = vc.genotypes
                if(gt.het) {
                    return new PopulationAlleleCounts(1,0,1)
                }
                else
                if(gt.homVar) {
                    return new PopulationAlleleCounts(sampleAN,1,0)
                }
                else 
                if(gt.homRef) {
                    return PopulationAlleleCounts.ZERO_COUNTS
                }
                else {
//                    assert false : "Unexpected / invalid genotype at position $vc.start ($vc)"
                    // This occurs if the genotype is './.' which means unascertained
                    return PopulationAlleleCounts.ZERO_COUNTS
                }
            }.sum()
        } 
    }

    /**
     * Generate a header for the population VCF. This is accomplished by copying headers from the first VCF and
     * adding INFO declarations for the population statistics entries.
     */
    void printVCFHeader() {
        VCF result = new VCF()
        result.headerLines = new VCF(opts.arguments()[0]).headerLines
        result.addInfoHeaders([
            '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">',
            '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">',
            '##INFO=<ID=GTC,Number=G,Type=Integer,Description="GenoType Counts. For each ALT allele in the same order as listed = 0/0,0/1,1/1,0/2,1/2,2/2,0/3,1/3,2/3,3/3">'
        ])
        result.headerLines[-1] = result.headerLines[-1].tokenize('\t')[0..7].join('\t')
        result.print(out)
    }

    /**
     * Infer the sexes for the input VCFs
     */
    void inferSexes() {
        
        Map<String,gngs.Sex> providedSexes = [:]
        if(opts.sexs)
             providedSexes = opts.sexs.collectEntries { def parts = it.tokenize(':'); [parts[0], gngs.Sex.decode(parts[1])] }
        
        List<VCF> vcfs = opts.arguments().collect { new VCF(it) }
        
        log.info "Computing sexes for ${vcfs*.samples*.size().sum()} samples ..."
        sampleSexes = vcfs.collect { VCF vcf ->
            vcf.samples.collect { String sampleId ->
                if(sampleId in providedSexes)
                    [sampleId, providedSexes[sampleId] ]
                else {
                    try {
                        [sampleId, vcf.guessSex(vcf.samples.indexOf(sampleId))]
                    }
                    catch(Exception e) {
                        throw new RuntimeException("Unable to infer sex for sample ${sampleId}. Please provide sex for this sample manually using -sex", e)
                    }
                }
            }
        }.sum().collectEntries()
        sexes = sampleSexes*.value
        log.info "Sample sexes: $sampleSexes"
    }
    
    static void main(String [] args) {
        cli('CreatePopulationStatisticsVCF <vcf1> <vcf2> <vcf3> ...', args) {
            o 'Output file', args:1, required: false
            sex 'Provide sex for a sample in the form <sampleid>:<MALE|FEMALE>', args: Cli.UNLIMITED
        }
    }

}

@CompileStatic
class PopulationAlleleCounts {
    int ac
    int hom
    int het
    
    final static ZERO_COUNTS = new PopulationAlleleCounts(0,0,0)
    
    PopulationAlleleCounts(final int an, final int hom, final int het) {
        this.ac = an; this.hom = hom; this.het = het;
    }
    
    PopulationAlleleCounts plus(PopulationAlleleCounts other) {
        if(other.is(ZERO_COUNTS))
            return this
        new PopulationAlleleCounts(ac+other.ac, hom+other.hom, het+other.het)
    }
}
