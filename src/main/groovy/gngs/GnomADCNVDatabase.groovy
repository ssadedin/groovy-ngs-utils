/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) Simon Sadedin, ssadedin<at>gmail.com
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

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import htsjdk.samtools.util.IntervalTree

/**
 * A {@link CNVDatabase} implementation for CNVs from the gnomAD
 * sv dataset.
 * <p>
 * Frequency statistics and overlaps are include all SVs (not just CNVs)
 * except unresolved breakends and variants failing the PASS filter.
 * 
 * @author simon.sadedin@gmail.com
 */
@CompileStatic
class GnomADCNVDatabase extends CNVDatabase {
    
    VCFIndex gnomad
    
    double mutualOverlapFraction = 0.5d

    Map<String, IntervalTree<Variant>> contigTrees = [:]
        
    public GnomADCNVDatabase(String path) {
        super();
        
        this.gnomad = new VCFIndex(path)

        VCF.parse(path) { Variant v ->
            IntervalTree<Variant> tree = contigTrees[v.chr]
            if(tree == null) {
                tree = contigTrees[v.chr] = new IntervalTree()
            }
            tree.put(v.pos, v.pos + v.size(), v)
            false
        }
    }
 
    public GnomADCNVDatabase(VCFIndex gnomad) {
        super();
        this.gnomad = gnomad;
    }
    

    @Override
    public double maxFreq(Map options, final IRegion region) {
        String contigNoChr = convertChr(region)
        
        Region rRegion = new Region(region)
        
        Variant maxFreqVariant = 
            (Variant)gnomad.iterator(contigNoChr, region.range.from, region.range.to)
                  .grep { Variant v -> included(v) && rRegion.mutualOverlap(v) > mutualOverlapFraction }
                  .max { v -> ((String)((Variant)v).info.POPMAX_AF).toDouble() }
        
        if(maxFreqVariant == null)          
            return 0.0d
            
        return ((String)maxFreqVariant.info.POPMAX_AF).toDouble();
    }

    @Override
    public List<Region> queryOverlapping(IRegion r) { 
        final String contigNorm = convertChr(r)

        final List<Region> variants = contigTrees[r.chr]?.overlappers((int)r.range.from, (int)r.range.to).grep { IntervalTree.Node<Variant> node ->
            included((Variant)node.value)
        }
        .collect { v ->
            createRegion(((IntervalTree.Node<Variant>)v).value)
        }
        return variants
    }
    
    private final Region createRegion(final Variant v) {
        
       final int an = ((String)v.info.AN).toInteger()
       final int ac = ((String)v.info.AC).toInteger()
       
       final boolean isDup = isDup(v)
       final boolean isDel = isDel(v)
       
       assert isDup || isDel
        
       new Region(v.chr, v.pos, v.pos + v.size(), 
//           variant: v,  
           sampleSize:an,
           duplication_frequency: isDup ?  ac/an : 0,
           deletion_frequency: isDel ?  ac/an : 0,
           observedGains: isDup ?  ac : 0,
           observedLosses: isDel ?  ac : 0,
           qual: v.qual,
           varType: isDup ? 'Gain' : 'Loss' // only handle gains and losses for now
       )
    }

    private final boolean isDel(final Variant v) {
        return (v.alt == '<DEL>' || v.alt == '<CN0>')
    }

    private final boolean isDup(final Variant v) {
        return (v.alt == '<DUP>' || v.alt == '<INS>')
    }

    @Override
    public double maxFreq(final IRegion region) {
        return maxFreq(null, region)
    }

    private final String convertChr(final IRegion region) {
        return region.chr.startsWith('chr') ? region.chr : 'chr' + region.chr
    }
    
    private final boolean included(final Variant v) {
        return (v.alt != '<BND>') && (v.filter == 'PASS') && (isDup(v) || isDel(v))         
    }
}
