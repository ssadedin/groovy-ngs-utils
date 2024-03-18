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
 * A {@link CNVDatabase} implementation for CNVs from Decipher 
 * Developmental Delay project.
 * <p>
 * Note: the database is not parsed by the constructor, you must
 * call the {@link #parse()} method yourself before using any
 * query methods.
 * <p>
 * See https://decipher.sanger.ac.uk/files/downloads/population_cnv.txt.gz
 * 
 * @author simon.sadedin@gmail.com
 */
class DecipherCNVs extends CNVDatabase {
     
    /**
     * Columns from schema of DGV in UCSC table
     */
    final static List<String> DECIPHER_COLUMNS = 
        ['population_cnv_id','chr','start','end','deletion_observations','deletion_frequency','deletion_standard_error','duplication_observations','duplication_frequency','duplication_standard_error','observations','frequency','standard_error','type','sample_size','study']
    
    String dddFile 
    
    @Delegate
    RangedData ddd
    
    DecipherCNVs(String dddFile) {
        this.dddFile = dddFile
    } 
    
    DecipherCNVs(RangedData ddd) {
        this.ddd = ddd
    } 
    
    DecipherCNVs parse() {
        this.ddd = new RangedData(dddFile, 1,2,3).load(readFirstLine: false, columnNames:DECIPHER_COLUMNS) { Region r ->
            r.setChr('chr' + r.chr)
            
            // Work around for seemingly incorrect parsing
            r.deletion_frequency = r.deletion_frequency.toDouble()
            r.duplication_frequency = r.duplication_frequency.toDouble()
            
            // These make it have properties compatible with DGV
            r.observedGains = r.duplication_observations
            r.observedLosses = r.deletion_observations
            r.sampleSize = r.sample_size
            if(r.type > 0)
                r.varType = 'Gain'
            else
            if(r.type < 0)
                r.varType = 'Loss'
            else
                r.varType = 'Gain+Loss'
        }
        return this
    }
    
    /**
     * Convenience method to return CNVS overlapping the specified region
     * 
     * @param r
     * @return
     */
    @Override
    List<Region> queryOverlapping(IRegion r) {
       return ddd.getOverlaps(r)*.extra
    }
    
    /**
     * Find the maximum frequency of this CNV within any study within DDD where the 
     * study has more than a minimum threshold size (default: 10 people).
     * 
     * @param region    region to search
     * @return  maximum frequency, or zero if no CNVs found
     */
    @Override
    double maxFreq(Map options=[:], IRegion region) {
        int minSampleSize = options.minSampleSize?:10
        List<Region> overlappingEntries = this.queryOverlapping(region)
        
        Region maxEntry = overlappingEntries.grep {
            (it.sample_size > minSampleSize) 
        }.max {
            return (it.deletion_frequency + it.duplication_frequency) 
        }
        
        if(maxEntry == null)
            return 0.0d
            
        return (maxEntry.deletion_frequency + maxEntry.duplication_frequency)
    } 
}
