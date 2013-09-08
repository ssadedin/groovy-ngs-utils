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

import groovy.transform.CompileStatic;

/**
 * Indexes a set of ranges so that the ranges crossing any point
 * can be identified quickly and easily.
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class RangeIndex {
    
    /**
     * The actual index - a map from contig / chromosome to the breakpoints
     * formed by the ranges added to the index, with the list of ranges
     * valid for each interval after the breakpoint.
     */
    TreeMap<Integer,List<IntRange>> ranges = new TreeMap()
    
    void add(int startPosition, int endPosition, Object extra = null) {
        add(extra != null ? new GRange(startPosition, endPosition-1, extra) : new IntRange(startPosition, endPosition-1))
    }
        
//    @CompileStatic
    void add(IntRange newRange) {
        
        // Any existing range over start position?
        int startPosition = (int)newRange.from
        int endPosition = (int)newRange.to
        Map.Entry<Integer,List<IntRange>> lowerEntry = ranges.lowerEntry(startPosition+1)
        if(lowerEntry == null || lowerEntry.value.isEmpty()) { // No range over this position: create one
            ranges.put(startPosition,[newRange])
            
            // If there are any ranges that start before the end of this new range,
            // add a breakpoint, and add this range to the intervening ranges
            int entryPos = startPosition
            Map.Entry higherEntry = ranges.higherEntry(startPosition)
            Map.Entry lastEntry = null
            while(higherEntry && higherEntry.key < endPosition) {
                higherEntry.value.add(newRange)
                lastEntry = higherEntry
                higherEntry = ranges.higherEntry(higherEntry.key)
            }
            
            // The last region needs to be split into two
            if(lastEntry) {
                if(!ranges.containsKey(endPosition)) {
                  List newRanges = new ArrayList(lastEntry.value.size())
                  newRanges.addAll(lastEntry.value)
                  ranges[endPosition] = newRanges
                }
            }
            else {
                // Need to mark the ending of the range so that future adds will 
                // look up the right point for adding new overlapping ranges
                ranges[endPosition] = []
            }
        }
        else { // Already a range over this position: we have to split it
            
            // TODO: does this have to be "max" ... or not?
            if(lowerEntry.value == null) {
                println "NULL VALUE"
            }
            else
            if(lowerEntry.value[0] == null) {
                println "VALUE HAS NULL CONTENTS"
            }
            
            int lowerEndPosition =  (int)lowerEntry.value[0].to
            
             // Split the beginning overlapping range
            // Copy the existing region to get all it's existing covered areas
            if(!ranges.containsKey(startPosition)) {
                List<IntRange> lowerSplitRegion = []
                if(lowerEndPosition > startPosition)
                    lowerSplitRegion.addAll(lowerEntry.value)
//                else 
//                    ranges[lowerEndPosition] = []
               
                // Add our new range as covered by the split part
                lowerSplitRegion.add(newRange)
                ranges[startPosition] = lowerSplitRegion
            }
            else {
                ranges[startPosition].add(newRange)
            }
                
            Map.Entry containedEntry = ranges.higherEntry(startPosition)
            if(!containedEntry) {
                // Can only happen when the lower region contains a single or 
                // set of identical regions
                ranges[lowerEndPosition] = [newRange]
            }
            
            List<Integer> rangesToAddTo = []
            while(containedEntry && containedEntry.key < endPosition) {
                
                Map.Entry higherEntry = ranges.higherEntry(containedEntry.key)   
                
                // Note: boundaryPoint is -1 when the ranged overlapped is an "ending" range - one at the 
                // border of a gap with no overlapping ranges. In that case we need to add ourselves to the 
                // boundary breakpoint and break
                int boundaryPoint = (higherEntry!=null) ? higherEntry.key : ((Integer)containedEntry.value[0]?.to?:-1)
                
                // If existing range is entirely contained within the one we are adding
                // then just add our new range to the list of ranges covered
                if(endPosition > boundaryPoint) { // Entirely contained
                    
                    // If there's no higher entry at all, then we can just add a new boundary
                    // and a region with only our range and break
                    if(higherEntry == null && boundaryPoint>=0) {
                        ranges[boundaryPoint] = [newRange]
                    }
                }
                else { // The start is contained, but the end is not : split the upper range
                    // Make a new range from the end of our range to the start of the next higher range
                    // It needs to have the previous range only, not our new range
                    // NOTE: list.clone() causes static compilation to fail with verify error
                    List<IntRange> clonedList = new ArrayList<IntRange>()
                    clonedList.addAll(containedEntry.value)
                    ranges[endPosition] = clonedList
                }
                rangesToAddTo << containedEntry.key
                
                containedEntry = higherEntry
            }
            rangesToAddTo.each { int startPos -> ranges[startPos] << newRange }
        }
    }
    
    List<Range> startingAt(int pos) {
        
        List<IntRange> result = this.ranges[pos]
        if(!result)
            return []
        
        // The given ranges could be starting at OR ending at the given position
        result.grep { it.from == pos }
    }
    
    List<Range> endingAt(int pos) {
        Map.Entry entry = this.ranges.lowerEntry(pos)
        if(!entry)
            return []
        List<IntRange> result = entry.value
        if(!result)
            return []
         result.grep { it.to == pos }
    }
    
    /**
     * Support for 'in' operator
     */
    boolean isCase(int position) {
        Map.Entry lowerEntry = ranges.lowerEntry(position+1)
        if(!lowerEntry)
            return false
        return lowerEntry.value.any { it.containsWithinBounds(position) }    
    }
    
    List<Range> getOverlaps(int position) {
        Map.Entry lowerEntry = ranges.lowerEntry(position+1)
        if(!lowerEntry)
            return []
        return lowerEntry.value.grep { it.containsWithinBounds(position) }
    }
    
    List<Range> intersect(int start, int end) {
        IntRange interval = start..end-1
        List<Range> result = []
        Map.Entry entry = ranges.lowerEntry(start+1)
        if(!entry) 
            entry = ranges.higherEntry(start)
           
        // Iterate as long as we are in the range
        while(entry != null && entry.key <= end) {
            for(List ranges in entry.value) {
               result << ranges 
            }
            entry = ranges.higherEntry(entry.key)
        }
        
        return result.collect { Math.max(it.from, start)..Math.min(it.to, end)}
    }    
    
    List<Range> subtractFrom(int start, int end) {

        List<Range> result = []
        List<Range> cutOut = intersect(start,end)       
        
        Range lastRange = start..start
        for(Range r in cutOut) {
            if(lastRange.to < r.from) {
              result.add(lastRange.to..r.from-1)
            }
            lastRange = r
        }
        if(lastRange.to<end)
          result.add(lastRange.to..(end-1))
        return result
    }    
    
    void dump() {
        // Get the starting point of all the regions
        def regions = ranges.keySet() as List
        regions += ranges.lastEntry().value[0].to
        
        def sizes = []
        for(int i=0; i<regions.size(); ++i) {
            if(i>0) {
                int regionSize = (regions[i] - regions[i-1])
                sizes << regionSize
            }
        }
        
        // Maximum region size
        int maxSize = regions.max()
        
        // Scale down if necessary
        int width = sizes.sum()
        def printSizes = sizes
        println "Total width = $width"
        println "Sizes = $sizes"
        if(width > 80) {
            printSizes = sizes*.multiply( 80.0f/width).collect { Math.round(it).toInteger() }
        }
        
        // now print out a representation of the regions
        println("|"+[sizes,printSizes].transpose().collect { it[0].toString().center(it[1],"-") }.join("|") +"|")
//        println("|"+sizes.collect { "-" * it }.join("|")+"|")
    }
    
}