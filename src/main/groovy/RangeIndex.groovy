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

import java.util.Iterator;

import groovy.lang.IntRange;
import groovy.transform.CompileStatic;

/**
 * Indexes a set of ranges so that the ranges crossing any point
 * can be identified quickly and easily.
 * <p>
 * The RangeIndex models ranges by their breakpoints. Each time a new range
 * is inserted, the breakpoints are stored in a tree and the list of overlapping ranges
 * at each breakpoint is tracked for each breakpoint entry. This makes it 
 * easy to find overlapping ranges at any point. This class is the primary driver
 * for the {@link BED} class which adds support for contigs so that a whole
 * reference sequence containing multiple contigs (eg: chromosomes) can be
 * indexed.
 * <p>
 * The RangeIndex is built upon Groovy's built in {@link IntRange} class. An 
 * important aspect of this class is that by default its end point is included
 * in the range, while by default BED ranges are exclusive of the end point.
 * Thus care needs to be taken in translating between the two. When you add
 * a range to the index it is treated as a BED range: that is, the end point
 * is exclusive. When you retrieve it, the end position will be <i>inclusive</i>
 * ie: it will be decremented by one from what you inserted.
 * <p>
 * RangeIndex does not de-duplicate identical ranges. If you put in multiple 
 * identical ranges they will be separate stored and separately returned as outputs
 * for queries and iteration over the index.
 * <p>
 * <b>Usage</b>
 * <p>
 * Adding ranges is done via the <code>add</code> method:
 * <pre>
 * RangeIndex index = new RangeIndex()
 * index.add(5..10)
 * index.add(20..30)
 * index.add(25..35)
 * </pre>
 * The {@link RangeIndex} class implements the {@link Iterable} interface, 
 * allowing generic Groovy collections methods to work:
 * <code>
 * index.each { println "A range from $it.from-$it.to" }
 * assert index.grep { it.from > 15 } == 2
 * </code>
 * <p>
 * Design notes: RangeIndex stores intervals in TreeMap (a balanced tree), indexed
 * by both start position and end position. Thus there is an entry in the index for
 * every "breakpoint" (start and end) of every range in the index. The value at each 
 * breakpoint is a list of the ranges that overlap the breakpoint. For example, if a range 
 * a..b is added (inclusive of both a and b), an entry at <code>a</code> is created and
 * the list at index <code>a</code> will contain the range <code>a..b</code>. Another
 * entry at <code>b+1</code> will be created, which does not contain <code>a</code>. 
 * (If there are no other overlapping ranges, the entry at <code>b+1</code> will be
 * empty).
 * <p>
 * To minimise memory use, RangeIndex accepts plain IntRange objects into the index. However
 * frequently it is desirable to associate a range to some more data. For this purpose,
 * the GRange class extends IntRange but contains an <code>extra</code> field that 
 * can be used to store arbitrary data in a range object.
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class RangeIndex implements Iterable<IntRange> {
    
    /**
     * The actual index - a map from position to the list of ranges
     * covering the region after the position. There is an entry in this
     * index for each position where a range starts or ends.
     */
    TreeMap<Integer,List<IntRange>> ranges = new TreeMap()
    
    void add(int startPosition, int endPosition, Object extra = null) {
        add(extra != null ? new GRange(startPosition, endPosition-1, extra) : new IntRange(startPosition, endPosition-1))
    }
        
    @CompileStatic
    void add(IntRange newRange) {
       
        // Any existing range over start position?
        int startPosition = (int)newRange.from
        int endPosition = (int)newRange.to
        Map.Entry<Integer,List<IntRange>> lowerEntry = ranges.lowerEntry(startPosition+1)
        
        // Inserting a range at the lowest end is relatively straightforward, so 
        // handle it separately
        if(lowerEntry == null || lowerEntry.value.isEmpty()) { 
            addLowestRange(newRange)
            return
        }
        
        // Already a range preceding this position: we to check for overlap and maybe split it
        if(lowerEntry.value == null) {
            println "NULL VALUE"
        }
        else
        if(lowerEntry.value[0].is(null)) {
            println "VALUE HAS NULL CONTENTS"
        }
            
        if(ranges.containsKey(startPosition)) { // If starts at exactly same position, share its entry
            ranges[startPosition].add(newRange)
            checkRanges(startPosition)
        }
        else { // starts at a new position, have to find overlapping ranges and split them
                
            // Add all the overlapping regions from the adjacent lower region to
            // our new breakpoint
            List<IntRange> lowerSplitRegion =  lowerEntry.value.grep { IntRange lwrRange -> lwrRange.to > startPosition }
                
            // Add our new range as covered by the split part
            lowerSplitRegion.add(newRange)
            ranges[startPosition] = lowerSplitRegion
            checkRanges(startPosition)
        }
                
        Map.Entry containedEntry = ranges.higherEntry(startPosition)
        List<Integer> rangesToAddTo = []
        
        boolean fullyContained = containedEntry && containedEntry.key >= endPosition
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
                if(higherEntry == null) { 
                    if(boundaryPoint>=0) {
                        ranges[boundaryPoint] = [newRange]
                        checkRanges(boundaryPoint)
                    }
                }
            }
            else { // The start is contained, but the end is not : split the upper range
                // Make a new range from the end of our range to the start of the next higher range
                // It needs to have the previous range only, not our new range
                // NOTE: list.clone() causes static compilation to fail with verify error
                List<IntRange> clonedList = containedEntry.value.grep { IntRange range -> endPosition < range.to }
                ranges[endPosition+1] = clonedList
                checkRanges(endPosition+1)
            }
            assert containedEntry.key < endPosition
            
            rangesToAddTo << containedEntry.key
                
            containedEntry = higherEntry
        }
        
        rangesToAddTo.each { int startPos -> 
            ranges[startPos] << newRange 
            checkRanges(startPos)
        }
        
        if(!ranges.containsKey(endPosition+1)) {
            if(fullyContained && lowerEntry) {
                ranges[endPosition+1] = lowerEntry.value.grep { IntRange range -> endPosition < range.to }
            }
            else {
                ranges[endPosition+1] = []
            }
        }
    }
    
//    @CompileStatic
    private void addLowestRange(IntRange newRange) {
        
        int startPosition = newRange.from
        int endPosition = newRange.to
        
        ranges.put(startPosition,[newRange])
//        checkRanges(startPosition)
        
        // If there are any ranges that start before the end of this new range,
        // add a breakpoint, and add this range to the intervening ranges
        int entryPos = startPosition
        Map.Entry<Integer, List<IntRange>> higherEntry = ranges.higherEntry(startPosition)
        Map.Entry<Integer, List<IntRange>> lastEntry = null
        while(higherEntry && higherEntry.key < endPosition) {
            higherEntry.value.add(newRange)
//            checkRanges(higherEntry.key)
            lastEntry = higherEntry
            higherEntry = ranges.higherEntry(higherEntry.key)
        }
        
        // The last region needs to be split into two
        if(lastEntry) {
            if(!ranges.containsKey(endPosition+1)) {
              List newRanges = lastEntry.value.grep { endPosition+1 in it }
              ranges[endPosition+1] = newRanges
//              checkRanges(endPosition+1)
            }
        }
        else {
            // Need to mark the ending of the range so that future adds will
            // look up the right point for adding new overlapping ranges
            if(!ranges.containsKey(endPosition+1)) {
                // println "Adding ${endPosition+1} to index ..."
                ranges[endPosition + 1] = []
            }
        }
    }
    
    @CompileStatic
    void checkRanges(int position) {
        return;
        
        List<IntRange> rangesAtPosition = ranges.get(position)
        
        
        if(rangesAtPosition == null)
            return 
           
        for(IntRange r in rangesAtPosition) { 
            assert r.from <= position && r.to>=position : "Range $r.from - $r.to is added incorrectly to position $position" 
        }
    }
    
    /**
     * Remove an existing range from the index. The existing range
     * MUST already be in the index, and cannot be a user-supplied range. 
     */
    void remove(Range r) {
        // There should be an entry where this range starts
        int lastPos = r.from-1
        List toRemove = []
        // Remove it from each breakpoint until we hit the end
        while(lastPos < r.to) {
          Map.Entry startEntry = this.ranges.higherEntry(lastPos)
          if(startEntry && startEntry.key > r.to)
              break
          
          if(!startEntry && lastPos == r.from-1)
              throw new IllegalArgumentException("Range $r.from-$r.to was not indexed, so cannot be removed")
          
          if(!startEntry)
              break
              
          List<IntRange> breakPointRanges = startEntry.value
//          println "Ranges are $breakPointRanges"
          // note issue when a new range starts at the same place as the range to be removed ends:
          // the new range probably is supposed to have the old range in its list? but it doesn't
          // so the last clause below stops us throwing on that case
          if(!breakPointRanges.isEmpty() && (!breakPointRanges.removeAll { it.from == r.from && it.to == r.to } && startEntry.key<r.to))
              throw new IllegalArgumentException("Range $r.from-$r.to is not present in the index, so cannot be removed")
          
          if(breakPointRanges.isEmpty()) {
              toRemove.add(startEntry.key)
          }
          lastPos = startEntry.key
        }
        
        toRemove.each { this.ranges.remove(it) }
        if(this.ranges.containsKey(r.to+1) && this.ranges[r.to+1].isEmpty()) {
            this.ranges.remove(r.to+1)
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
    
    /**
     * Return a list of ranges stored in the index that overlap the specified
     * range. Both ends of the range are *inclusive*.
     * 
     * @param start first position to look for overlaps
     * @param end   last position to look for overlaps (inclusive)
     * @return  List of ranges overlapping the specified start -> end range
     */
    List<Range> getOverlaps(int start, int end) {
        getOverlaps(start,end,false)
    }
    
    /**
     * Return true if the given range overlaps at least one range in this
     * RangeIndex.
     * 
     * @param start start of range to test (inclusive)
     * @param end   end point of range to test (inclusive)
     * @return  true iff at least one range in this index overlaps the given start and end point
     */
    boolean overlaps(int start, int end) {
        !getOverlaps(start,end, true).isEmpty()
    }
    
    @CompileStatic
    List<Range> getOverlaps(int start, int end, boolean returnFirst) {
        IntRange interval = start..end-1
        List<Range> result = []
        Map.Entry<Integer,List<IntRange>> entry = ranges.lowerEntry(start+1)
        if(!entry) 
            entry = ranges.higherEntry(start)
            
        // Iterate as long as we are in the range
        while(entry != null && entry.key <= end) {
            for(IntRange r in entry.value) {
               if(r.from<=end && r.to>=start) {
                   result.add(r)
                   if(returnFirst)
                       return result
               }
            }
            entry = ranges.higherEntry(entry.key)
        }
        return result
    }
    
    /**
     * Return a list of ranges that intersect the given start and end points.
     * <p>
     * <em>Note:</em>The start and end point are both treated as <b>inclusive</b>.
     * @param start start of range to find intersections for
     * @param end   end of range to find intersections for
     * @return  List of ranges that intersect the specified range.
     */
    List<Range> intersect(int start, int end) {
       def result = getOverlaps(start,end)
       return result.collect { Math.max(it.from, start)..Math.min(it.to, end)}
    }    
    
    /**
     * Subtract all the ranges in this range index from the given range and return
     * the resulting list of ranges.
     * <p>
     * <em>Note:</em>both start and end are considered <b>inclusive</b>.
     * 
     * @param start
     * @param end
     * @return  List of ranges left after the ranges in this RangeIndex are removed from 
     *          the specified region.
     */
    List<Range> subtractFrom(int start, int end) {

        List<Range> result = []
        List<Range> cutOut = intersect(start,end)       
        
        // None of the regions in our bed file overlap the range specified,
        // so just return the whole range
        if(cutOut.isEmpty())
            return [start..end-1]
        
        Range lastRange 
        if(start < cutOut[0].from)
            lastRange = start-1..start-1
        else {
            lastRange = cutOut[0]
        }
            
        for(Range r in cutOut) {
            if(lastRange.is(r))
                continue
            if(lastRange.to < r.from) {
              result.add((lastRange.to+1)..r.from-1)
            }
            lastRange = r
        }
        if(lastRange.to<end) {
          // Important note: Groovy IntRange will reverse to and from
          // if from > to. So, we can't construct the range until after we 
          // do the test here about whether to > from
          // alternatively we could have used rangeToAdd.@to / rangeToAdd.@from
          int fromPos = (lastRange.to+1)
          int toPos = (end-1)
          if(toPos >= fromPos) {
              IntRange rangeToAdd = fromPos..toPos 
              result.add(rangeToAdd)
          }
        }
        return result
    }    
    
    Range nextRange(int pos) {
        def entry = [key: pos]
        while(true) {
          entry = ranges.higherEntry(entry.key)
          if(!entry)
              return null
          List<IntRange> nextRanges = entry.value.grep { it.from >= pos } 
          if(nextRanges)
              return nextRanges[0]
        }
    }
    
    Range previousRange(int pos) {
        def entry = [key: pos]
        while(true) {
          entry = ranges.lowerEntry(entry.key)
          if(!entry)
              return null
          List<IntRange> nextRanges = entry.value.grep { it.to <= pos } 
          if(nextRanges)
              return nextRanges[0]
        }
    }
    
    Range nearest(int pos) {
        List<IntRange> overlaps = this.getOverlaps(pos)
        if(overlaps) {
            return overlaps.min {Math.min(it.to-pos, pos-it.from)}
        }
        else {
            IntRange prv = previousRange(pos)
            IntRange nxt = nextRange(pos)
            if(!prv)
                return nxt
            if(!nxt)
                return prv
            
            return (pos-prv.to) < (nxt.from - pos) ? prv : nxt
        }
    }
    
    Range first() {
        if(ranges.isEmpty())
            return null
            
        ranges.firstEntry().value[0]
    }
    
    Range last() {
        if(ranges.isEmpty())
            return null
        Map.Entry last = ranges.lastEntry()
        while(!last.value)
            last = ranges.lowerEntry(last.key)
            
        return last.value[0]    
    }
    
    void dump() {
        // Get the starting point of all the regions
        def regions = ranges.keySet() as List
        def lastRange = ranges.lastEntry().value[0]
        if(lastRange)
          regions += lastRange.to
        
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
        int width = sizes?sizes.sum():0
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

    @Override
    public Iterator<IntRange> iterator() {
        iteratorAt(-1)
    }
    
    /**
     * Return an iterator that will return each range in the index
     * in genomic order by starting position, starting from the given 
     * position.
     */
    public Iterator<IntRange> iteratorAt(int startingPos) {
        
        if(startingPos >= 0) {
            startingPos = ranges.containsKey(startingPos) ? startingPos : ranges.higherKey(startingPos)
        }
		int startingPosTmp = startingPos
        return new Iterator<IntRange>() {
            Integer pos = startingPosTmp
            int index = 0
            List<IntRange> activeRanges = []
            IntRange nextRange = null
            
            boolean hasNext() {
                if(nextRange)
                    return true
                    
                // Note: pos becomes null when iteration through index exhausted
                // nextRange is set null at each call of next()
                while(pos != null && (nextRange==null || nextRange.from < pos)) {
                  nextRange = findNext()
                }
                return nextRange != null
            }
            
            IntRange next() {
                
                if(!hasNext()) // Note: populates nextRange
                    throw new IndexOutOfBoundsException("No more ranges in this index")
                
                IntRange result = nextRange
                nextRange = null
                return result
            }
            
            IntRange findNext() {
                if(pos == -1) {
                    pos = ranges.firstKey()
                }
                else
                if(index > ranges[pos].size()-1) {
                    index = 0
                    
                    pos = ranges.higherKey(pos)
                    while(pos != null && ranges[pos].isEmpty()) {
                      pos = ranges.higherKey(pos)
                    }
                }
                
                if(pos == null) 
                    return null
                
                return ranges[pos][index++]
            }
            
            void remove() {
                throw new UnsupportedOperationException()
            }
        }
    }
    
    public Iterator<IntRange> reverseIterator() {
        reverseIteratorAt(-1)
    }
    
    /**
     * Return an iterator that will return each range that starts
     * at or before the given position, proceeding backwards
     * through the genome, ordered by start position.
     * 
     * @author simon.sadedin@mcri.edu.au
     */
    public Iterator<IntRange> reverseIteratorAt(int startingPos) {
        
		int startingPosTmp = startingPos
        if(startingPos >= 0) {
            // We have to ensure to start at least 1 higher than
            // the first entry we want to return, because the findNext()
            // will call index.lowerKey()
            if(ranges.containsKey(startingPos)) {
                startingPosTmp = startingPos+1
            }
        }
        
        return new Iterator<IntRange>() {
            
            // Genomic position in range index
            Integer pos = startingPosTmp
            
            // position in values at current range index
            int index = -1
            
            List<IntRange> activeRanges = []
            IntRange prevRange = null
            
            boolean hasNext() {
                if(prevRange)
                    return true
                    
                // Note: pos becomes null when iteration through index exhausted
                // prevRange is set null at each call of next()
                while(pos != null && (prevRange==null || prevRange.from < pos)) {
                  
                  if(prevRange!=null && prevRange.from < pos) {
                      // println "Skip $prevRange.from - $prevRange.to because ${prevRange.from} < $pos"
                  }
                  prevRange = findNext()
                  
                }
                return prevRange != null
            }
            
            
            IntRange next() {
                
                if(!hasNext()) // Note: populates prevRange
                    throw new IndexOutOfBoundsException("No more ranges in this index")
                
                IntRange result = prevRange
                prevRange = null
                return result
            }
            
            IntRange findNext() {
                if(pos == -1) {
                    pos = ranges.lastKey()
                    if(pos != null) {
                        index = ranges[pos].size()-1
                        if(index<0)   // This actually forces the return below to return null (in groovy, [0] == null)
                            index = 0 // which forces the iteration to occur again, to find the right key
                    }
                }
                else
                if(index < 0) {
                    pos = ranges.lowerKey(pos)
                    while(pos != null && ranges[pos].isEmpty()) {
                      pos = ranges.lowerKey(pos)
                    }
                    
                    if(pos != null) {
                        index = ranges[pos].size()-1
                    }
                }
                
                if(pos == null)
                    return null
                
                return ranges[pos][index--]
            }
            
            void remove() {
                throw new UnsupportedOperationException()
            }
        }
    }
    
    /**
     * Merge all overlapping ranges together to make simplified regions
     * representing all the ranges covered by any range in this RangeInde.
     * 
     * @return
     */
    RangeIndex reduce() {
        RangeIndex reduced = new RangeIndex()

        // We take advantage of the fact that we iterate the ranges in order of genomic start position
        IntRange currentRange = null
        for(IntRange r in this) {
            if(currentRange == null) {
                currentRange = r
                continue
            }
            
            assert r.from >= currentRange.from
            
            if(GRange.overlaps(r,currentRange)) { 
                currentRange = Math.min(r.from, currentRange.from)..Math.max(r.to, currentRange.to)
            }
            else {
                reduced.add(currentRange)
                currentRange = r
            }
        }
        if(currentRange != null)
            reduced.add(currentRange)
            
        return reduced
    }
    
    /**
     * Returns a range for each boundary within the index.
     * 
     * @param c
     * @return
     */
    @CompileStatic
    Iterator<GRange> getBoundaries() {

        Iterator<Map.Entry<Integer,List<IntRange>>> i = ranges.iterator()

        return new Iterator<GRange>() {
            
            Map.Entry<Integer,IntRange> prev = null
            
            IntRange prevRange = null
            
            boolean hasNext() {
//                println "hasNext = " + i.hasNext() + " prev = " + prev
                i.hasNext() || (prev != null && !prev.value.empty)
            }
            
            GRange next() {
                
                Map.Entry<Integer,IntRange> curr = null
                if(prev == null) {
                    prev = i.next()
                }
                
                curr = i.hasNext() ? i.next() : prev
                
                
                // When we reach a gap we skip to the next range
                // as we do not output ranges covering gaps
                // A gap is indicated by a previous breakpoint where
                // the list of overlapping ranges is empty
                if(prev.value.empty && i.hasNext()) {
                    prev = curr
                    prevRange = null
                    return next()
                }
                
                GRange range 
                if(curr != prev) {
                    
                    int endPos = curr.key
                    int startPos =  prev.key                    
                    int beforeEndPos = endPos - 1
                    

                    // We want to output contiguous ranges, ie: the end of
                    // of the previous range is the start of the next,
                    // EXCEPT if there is a gap. For these gaps, we only know
                    // for real that it is a gap from information available in the 
                    // previous iteration. So in the previous one, we set prevRange
                    // non-null, only if its 'to' should be used as a starting point
                    if(prevRange != null) {
                        startPos = prevRange.to
                    }
                    
//                    println "=" * 80
//                    println "PrevRange=${prevRange?.from}-${prevRange?.to}"
//                    println "Curr = $curr"
//                    println "Prev = $prev"
//                    
                    // Is this an ending position? If so, return the pevious value
                    if(prev.value.any { it.to == beforeEndPos }) {
                        endPos -= 1
                    }
                    
                    range = new GRange(startPos, endPos, curr.value)
                    
                    // Only store the current output as prevRange if the current range is not 
                    // the beginning of a gap. It's a gap if there are no overarching / overlapping
                    // ranges, which would be stored in its value. That is true if every range
                    // that is in the previous position ends at the position prior to our current position.
                    if(prev.value.every {it.to == beforeEndPos}) {
                        prevRange = null
                    }
                    else {
                        prevRange = range
                    }
                    
                    prev = curr
                }
                else {
                    range = new GRange(prev.key, curr.value.max { it.to }.to, curr.value)
                    prev = null
                }
                return range
            }
            
            void remove() {
                throw new UnsupportedOperationException()
            }
        }
 
    }
}