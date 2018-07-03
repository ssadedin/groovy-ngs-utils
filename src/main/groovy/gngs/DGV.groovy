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

import java.util.zip.GZIPInputStream

import groovy.transform.CompileStatic

/**
 * A {@link CNVDatabase} implementation for CNVs from the Database of
 * Genomic Variants (sourced via UCSC).
 * <p>
 * Note: the database is not parsed by the constructor, you must
 * call the {@link #parse()} method yourself before using any
 * query methods.
 * 
 * @author simon.sadedin@gmail.com
 */
class DGV extends CNVDatabase {
    
    /**
     * Columns from schema of DGV in UCSC table
     */
    final static List<String> DGV_COLUMNS = ["bin", "chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "varType", "reference", "pubMedId", "method", "platform", "mergedVariants", "supportingVariants", "sampleSize", "observedGains", "observedLosses", "cohortDescription", "genes", "samples"]
    
    String dgvFile 
    
    @Delegate
    RangedData dgv
    
    DGV(String dgvFile) {
        this.dgvFile = dgvFile
    } 
    
    DGV(RangedData dgv) {
        this.dgv = dgv
    } 
     
    /**
     * Parse the data source set for this database
     */
    DGV parse() {
        this.dgv = new RangedData(dgvFile, 1,2,3).load(columnNames:DGV_COLUMNS) 
        return this
    }
    
    /**
     * return CNVS overlapping the specified region
     * @param r
     * @return
     */
    List<Region> queryOverlapping(Region r) {
       return dgv.getOverlaps(r)*.extra
    }
    
    /**
     * Find the maximum frequency of this CNV within any study within DGV where the 
     * study has more than a minimum threshold size (default: 10 people).
     * 
     * @param region    region to search
     * @return  maximum frequency, or zero if no CNVs found
     */
    double maxFreq(Map options=[:], Region region) {
        int minSampleSize = options.minSampleSize?:10
        List<Region> overlappingEntries = this.queryOverlapping(region)
        
        
        Region maxEntry = overlappingEntries.grep {
            (it.sampleSize > minSampleSize) && (it.observedGains + it.observedLosses < it.sampleSize)
        }.max {
            (it.observedGains + it.observedLosses) / it.sampleSize
        }
        
        if(maxEntry == null)
            return 0.0d
            
        return (maxEntry.observedGains + maxEntry.observedLosses) / maxEntry.sampleSize
    } 
    
   List<Region> filterByType(String type, List<Region> regions) {
       if(type == "LOSS")
          return regions.grep { it.varType in ["Loss","Gain+Loss","loss","gain+loss","deletion"] }
       else
          return regions.grep { it.varType in ["Gain","Gain+Loss","gain","gain+loss","duplication"] }
   }
    
}
