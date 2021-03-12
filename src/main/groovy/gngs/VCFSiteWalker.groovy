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
package gngs

import java.util.List

import groovy.transform.CompileStatic
import groovy.transform.stc.ClosureParams
import groovy.transform.stc.FromString
import groovy.transform.stc.SimpleType
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader

/**
 * A datastructure to hold thhe state of a given VCF that is being walked
 * 
 * @author Simon Sadedin
 */
@CompileStatic
class VCFWalkPosition {
    
    Iterator<VariantContext> iter
    
    List<VariantContext> variant
    
    VariantContext pending = null
    
    void next() {
        VariantContext v = nextImpl()
        if(v.is(null)) {
            variant = null
            return
        }
        variant = [v]
        final long start = v.start
        while(pending?.start == v.start) {
            variant.add(nextImpl())
        }
    }
    
    /**
     * Lets us peek at the next variant coming up
     * @return  the next variant from the iterator or null if there is none
     */
    VariantContext nextImpl() {
        if(pending == null) {
            if(iter.hasNext())
                pending = iter.next()
        }
        
        if(pending == null)
            return null
        
        VariantContext result = pending
        if(iter.hasNext()) {
            pending = iter.next()        
        }
        else {
            pending = null
        }
        return result
    }
}

/**
 * Traverses a set of VCFs site by site, such that all sites in the VCFs are visited, or alternatively,
 * allele by allele (such at that each site is processed multiple times, once for each ref/alt allele combination).
 * <p>
 * <b>Note: </b> Requires sorted, normalised VCF
 * <p>
 * Example usage:
 * <pre>
 * walker = new VCFSiteWalker(["test1.vcf","test2.vcf","test3.vcf"])
 * walker.walk { List<VariantContext> variantsAtSite ->
 *    println "There are ${variantsAtSite.size()} variants at ${variantsAtSite[0].start}"
 * }
 * </pre>
 * 
 * @author Simon Sadedin
 */
class VCFSiteWalker {
    
    List<VCFWalkPosition> vcfs
    
    public VCFSiteWalker(List<String> vcfs) {
        super();
        this.vcfs = vcfs.collect { 
            new VCFWalkPosition(iter: new VCFFileReader(new File(it), false).iterator())
        }
    }    
    
    public VCFSiteWalker(List<String> vcfs, List<String> names) {
        super();
        this.vcfs = [vcfs,names].transpose().collect { vcf, name ->
            VCFFileReader reader = new VCFFileReader(new File(vcf), false)
            reader.reader.codec.name = name
            new VCFWalkPosition(iter: reader.iterator())
        }
    }    
 
    
    /**
     * The primary entry point for using the VCFSiteWalker. Call this method passing a closure 
     * as a callback to receive each VCF site as a list if {@link VariantContext} objects.
     * 
     * @param callback  Closure as callback
     */
    @CompileStatic
    void walk(@ClosureParams(value=FromString, options='java.util.List<htsjdk.variant.variantcontext.VariantContext>') Closure callback) {
        
        vcfs*.next()
        
        while(vcfs.any { it.variant != null }) {
            processNextPosition(callback)
        }
    }

    /**
     * Walk all the input VCFs by allele - that is, the callback is invoked once for each
     * unique position/allele within the VCF
     * <p>
     * Note: this differs from the basic {@link #walk} method because there can
     * be multiple alt alleles at the same site. In that method, all the alleles at the same
     * site will be included together in the same callback. 
     * <p>
     * There can also be multiple reference alleles (eg: indel from one individual overlapping
     * a SNP from a second). Therefore this method further groups the variants into distinct 
     * ref/alt groups, and each group is output separately as a separate VCF record.
     * 
     * @param callback  Closure as callback
     */
    @CompileStatic
    void walkByAllele(@ClosureParams(value=FromString, options='java.util.List<htsjdk.variant.variantcontext.VariantContext>') Closure callback) {
        
        vcfs*.next()
        
        while(vcfs.any { it.variant != null }) {
            processNextPosition { List<VariantContext> variants ->
                processByAllele(variants, callback)
            }
        }
    }
    
    @CompileStatic
    private void processByAllele(List<VariantContext> variants, @ClosureParams(value=FromString, options='java.util.List<htsjdk.variant.variantcontext.VariantContext>') Closure callback) {
        Map<String, List<VariantContext>> alleles = variants.groupBy { VariantContext ctx ->
            ctx.getReference().baseString + '/' + ctx.getAlternateAllele(0).baseString
        }

        alleles.each { String bases, List<VariantContext> alleleVariants ->
            callback(alleleVariants)
        }
    }
    
    /**
     * Process a single site and advance the VCF iterators used to the next position
     * 
     * @param callback  
     */
    private void processNextPosition(Closure callback) {
        
        // Group by position, then order by those
        List<List<VCFWalkPosition>> nextVariants = vcfs.groupBy { VCFWalkPosition vcf ->
            List<VariantContext> v = vcf?.variant
            if(!v)
                return Long.MAX_VALUE
                
            VariantContext v0 = v[0]
            return XPos.computePos(v0.contig, v0.start)
        }
        .sort {
            it.key
        }.collect { 
            it.value // there is a list of variants from each VCF at each position: join them together
        }
            
        processPosition(nextVariants[0], callback)
    }
    
    private void processPosition(List<VCFWalkPosition> variants, Closure callback) {
        callback(variants*.variant.sum())
        variants*.next()
    }
}
