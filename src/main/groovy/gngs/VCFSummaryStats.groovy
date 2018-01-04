/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2017 Simon Sadedin, ssadedin<at>gmail.com
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
 * A set of statistics summarizing how a VCF was filtered
 * 
 * @author Simon Sadedin
 */
class VCFSummaryStats {
    
    int total = 0
    
    int excludeByDiff = 0
    
    int excludeByCons = 0
    
    int excludeByMaf = 0
    
    int excludeByTarget = 0
    
    int excludeComplex = 0
    
    int excludeNotPresent = 0
    
    int excludeByPreFilter = 0
    
    int excludeByFilter = 0
    
    int excludeByMasked = 0
    
    int totalIncluded
    
    @Override
    String toString() {
        String.valueOf("total=$total, incl=$totalIncluded, excl: cons=$excludeByCons, maf=$excludeByMaf, complex=$excludeComplex, prefilter=$excludeByPreFilter, filter=$excludeByPreFilter, masked=$excludeByMasked, target=$excludeByTarget, diff=$excludeByDiff")
    }
}

