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
    
    VariantContext variant
    
    void next() {
        if(iter.hasNext())
            variant = iter.next()
        else
            variant = null
    }
}

/**
 * Traverses a set of VCFs site by site, such that all sites in the VCFs are visited
 * one at a time.
 * <p>
 * <b>Note: </b> Requires sorted, indexed VCF
 * <b>Note2:</b> Sites may be repeated if they are repeated in a given VCF (eg: a VCF where multiallelic sites have been decomposed into primitives)
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
            new VCFWalkPosition(iter: new VCFFileReader(new File(it)).iterator())
        }
    }    
    
    /**
     * The primary entry point for using the VCFSiteWalker. Call this method passing a closure 
     * as a callback to receive each VCF site as a list if {@link VariantContext} objects.
     * 
     * @param callback  Closure as callback
     */
    void walk(@ClosureParams(value=SimpleType, options='gngs.VCFWalkPosition') Closure callback) {
        
        vcfs*.next()
        
        while(vcfs.any { it.variant != null }) {
            processNextPosition(callback)
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
            VariantContext v = vcf?.variant
            if(!v)
                return Long.MAX_VALUE
                
            return XPos.computePos(v.contig, v.start)
        }
        .sort {
            it.key
        }.collect { 
            it.value
        }
            
        processPosition(nextVariants[0], callback)
    }
    
    private void processPosition(List<VCFWalkPosition> variants, Closure callback) {
        callback(variants*.variant)
        variants*.next()
    }
}
