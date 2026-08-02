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

/**
 * A tool to create a "population statistics VCF" that represents the frequency of
 * variants from a given set of samples by entries in the info field.
 * <p>
 * The following population metrics are calculated:
 * <li>AC   -   the number of alleles that were observed at the site, binned according to frequencies [0-0.25, 0.25-0.5, 0.5-0.75, 0.75-1]
 * <li>AN   -   the total number of observed alleles
 * </p>
 * <b>Limitations:</b>
 * The current implementation has a number of limitations:
 * 
 * <li>Only VCF information is considered, not gVCF, which is to say, an allele is counted as 
 *     unobserved at a site even if there was no coverage and therefore was not ascertained at 
 *     that site.
 * 
 * These limitations should be corrected to generate truly accurate counts.
 * 
 * @author Macgregor Todd & Simon Sadedin
 */
@Log
class CreateMitoPopulationStatisticsVCF extends ToolBase {

    int numVCFs

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

            List<VCF> vcfs = opts.arguments().collect { new VCF(it) } // will help us compute the number of samples

            int numVCFs = vcfs*.samples*.size().sum()

            log.info "Computing allele count bins for $numVCFs samples ..."

            printVCFHeader(numVCFs)

            out.flush()

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
        // 

        Map<String, List<VariantContext>> alleles = 
            variants.groupBy { VariantContext ctx -> 
                ctx.getAlternateAllele(0).baseString
            }
            
        if(alleles.size()>1) {
            log.info "Site ${variants[0].contig}:${variants[0].start} is multiallelic"
        }
        
        for(List<VariantContext> allele : alleles.values()) {
            VariantContext v0 = allele[0]
            List<Integer> ac = computeAlleleCount(allele)
            int an = (int)ac.sum()
            printVCFSite(v0, ac, an)
        }
    }

    /**
     * Print out a line of the output VCF based on the given computed statistics
     * 
     * @param v0
     * @param ac
     * @param an
     */
    @CompileStatic
    protected void printVCFSite(VariantContext v0, List ac, int an) {
        String ac_forprint = ac.toString()[1..-2].replaceAll(' ', '')
        out.println([
            v0.contig,
            v0.start,
            v0.getID(),
            v0.reference.baseString,
            v0.alternateAlleles[0],
            '.',
            '.',
            "AC=$ac_forprint;AN=$an"
        ].join('\t'))
    }
    
    /**
     * Compute the allele count at a site, based on the given list of variants at that site.
     * 
     * @return  the number of alleles observed, based on the het / hom status of the variants
     *          and the sex of the respective samples
     */
    @CompileStatic
    List<Integer> computeAlleleCount(List<VariantContext> variants) {
        Set<String> samples = new HashSet(numVCFs*2)
        List<LinkedHashMap> ac_raw = (List)variants.collect { VariantContext vc -> // count across VCFs
            vc.genotypes.collect { Genotype gt ->     // sum across samples within this vcf // *** typically no trio vcfs in mito, but this closure maintains compatibility

                LinkedHashMap<String,Integer> bins = [
                    b1:0,
                    b2:0,
                    b3:0,
                    b4:0
                ]

                // Make us robust to samples being provided twice
                if(samples.contains(gt.sampleName in samples))
                    return bins
                else
                    samples.add(gt.sampleName)

                if((gt.extendedAttributes.AF as Double) <= 0.25) {
                    //log.info "$gt.extendedAttributes.AF into bin 1"
                    bins.b1 += 1
                }
                else if((gt.extendedAttributes.AF as Double) <= 0.5) {
                    bins.b2 += 1
                }
                else if((gt.extendedAttributes.AF as Double) <= 0.75) {
                    bins.b3 += 1
                }
                else {
                    bins.b4 += 1
                }
                return bins

            }.collectMany{it.entrySet()}.groupBy{it.key}.collectEntries{[it.key, it.value.sum{it.value}]} // maintains trio compatibility
        }
		// sum across VCFs
        LinkedHashMap<String,Integer> ac_summed = (LinkedHashMap)ac_raw.collectMany{it.entrySet()}.groupBy{it.key}.collectEntries{[it.key, it.value.sum{it.value}]}
        List<Integer> ac = ac_summed.collect{it.value}
    }

    /**
     * Generate a header for the population VCF. This is accomplished by copying headers from the first VCF and
     * adding INFO declarations for the population statistics entries.
     */
    void printVCFHeader(int numVCFs) {
        VCF result = new VCF()
        result.headerLines = new VCF(opts.arguments()[0]).headerLines
        result.addInfoHeaders([
            '##INFO=<ID=AC,Number=A,Type=List<Integer>,Description="Allele counts in genotypes, for each ALT allele, PER allele frequency bin. Frequency bins are [0-0.25, 0.25-0.5, 0.5-0.75, 0.75-1]">',
            '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">',
            "##These mito obs counts were calculated from $numVCFs VCFs"
        ])
        result.headerLines[-1] = result.headerLines[-1].tokenize('\t')[0..7].join('\t')
        result.print(out)
    }

    static void main(String [] args) {
        cli('CreateMitoPopulationStatisticsVCF <vcf1> <vcf2> <vcf3> ...', args) {
            o 'Output file', args:1, required: false
            args: Cli.UNLIMITED
        }
    }
}
