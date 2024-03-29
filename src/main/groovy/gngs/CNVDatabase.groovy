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

/**
 * Abstract interface that can be implemented by different
 * CNV sources.
 * <p>
 * CNVs are operated on as {@link Region} objects, but 
 * implementations will set implementation-specific properties
 * that can be accessed as properties directly on the class.
 * 
 * @author Simon Sadedin
 */
abstract class CNVDatabase {
    
    /**
     * Query for regions overlapping region {@link #r}
     * 
     * The returned regions are expected to have a number of properties added as expando properties:
     * 
     *   <li>deletion_frequency
     *   <li>duplication_frequency
     *   <li>observedGains
     *   <li>observedLosses 
     *   <li>sampleSize
     *   <li>varType  - one of ["Loss","Gain+Loss","loss","gain+loss","deletion","Gain","gain","duplication"]
     * 
     * @return list of CNVS overlapping the specified region
     */ 
    abstract List<Region> queryOverlapping(IRegion r) 
    
    /**
     * Find the maximum frequency of this CNV within any study within the database
     * where the study has more than a minimum threshold size (default: 10 people).
     * 
     * @param   region  region to search
     * @return  maximum frequency, or zero if no CNVs found
     */
    abstract double maxFreq(Map options=[:], IRegion region) 
}
