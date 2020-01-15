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

package gngs

import java.util.Iterator;

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
 * Design notes: RangeIndex stores intervals as a TreeMap (a balanced tree), indexed
 * by breakpoints. Thus there is an entry in the index for every "breakpoint" (start 
 * and end) of every range in the index. The value at each breakpoint is a list of the 
 * ranges that overlap the breakpoint. For example, if a range 
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
    
    private int numRanges = 0
    
    @CompileStatic
    public int getNumRanges() {
        return this.@numRanges
    }
    
    /**
     * Add the given range to this index
     * <p>
     * If the 'extra' argument is specified then a GRange will be added with the <code>extra</code>
     * as its extra data argument. Otherwise, a plain IntRange will be added.
     */
    @CompileStatic
    void add(int startPosition, int endPosition, Object extra = null) {
        add(extra != null ? new GRange(startPosition, endPosition-1, extra) : new IntRange(startPosition, endPosition-1))
    }
        
    @CompileStatic
    void add(IntRange newRange) {
        
        ++numRanges
       
        // Any existing range over start position?
        final int startPosition = (int)newRange.from
        final int endPosition = (int)newRange.to
        Map.Entry<Integer,List<IntRange>> lowerEntry = ranges.lowerEntry(startPosition+1)
        
        // Inserting a range at the lowest end is relatively straightforward, so 
        // handle it separately
        if(lowerEntry == null || lowerEntry.value.isEmpty()) { 
            addLowestRange(newRange)
            return
        }
        
        assert lowerEntry.value != null
        assert !lowerEntry.value[0].is(null) 
            
        // Already a range preceding this position: we to check for overlap and maybe split it
        if(ranges.containsKey(startPosition)) { // If starts at exactly same position, share its entry
            ranges[startPosition].add(newRange)
            checkRanges(startPosition)
        }
        else { // starts at a new position, have to find overlapping ranges and split them
            splitStartRange(newRange, lowerEntry)
        }
                
        Map.Entry containedEntry = ranges.higherEntry(startPosition)
        while(containedEntry && containedEntry.key < endPosition) {
            Map.Entry higherEntry = ranges.higherEntry(containedEntry.key)   
            addContainedEntry(containedEntry, higherEntry, newRange)
            containedEntry = higherEntry
        }
       
        if(!ranges.containsKey(endPosition+1)) {
            boolean fullyContained = containedEntry && containedEntry.key >= endPosition
            if(fullyContained && lowerEntry) {
                ranges[endPosition+1] = lowerEntry.value.grep { IntRange range -> endPosition < range.to }
            }
            else {
                ranges[endPosition+1] = []
            }
        }
    }
    
    /**
     * Add a reference to a range to a breakpoint that falls in the middle of the range
     * 
     * @param containedEntry    the entry for the breakpoint that is contained by the range
     * @param higherEntry       the entry for the next breakpoint. <code>null</code> if there is no higher entry.
     * @param newRange          the range to add to the breakpoint entry
     */
    @CompileStatic
    void addContainedEntry(final Map.Entry<Integer,List<IntRange>> containedEntry, final Map.Entry<Integer,List<IntRange>> higherEntry, final IntRange newRange) {
       
        final int endPosition = newRange.to
        
        final IntRange firstContainedRange = containedEntry.value[0]
        
        // Note: boundaryPoint is -1 when the ranged overlapped is an "ending" range - one at the
        // border of a gap with no overlapping ranges. In that case we need to add ourselves to the
        // boundary breakpoint and break
        // TODO: why is only the first range considered for this test?
        final int boundaryPoint = (higherEntry!=null) ? higherEntry.key : (firstContainedRange?.to?:-1)

        // If existing range is entirely contained within the one we are adding
        // then just add our new range to the list of ranges covered
        // 
        //   |-------------------------|  # Range to add
        //           |------|             # Other range whose boundary is inside this one
        if(endPosition > boundaryPoint) { // Entirely contained

            // If there's no higher entry at all, then we can just add a new boundary
            // and a region with only our range and break
            if(higherEntry == null) {
                if(boundaryPoint>=0) {
                    ranges[boundaryPoint] = [newRange]
                    // checkRanges(boundaryPoint)
                }
            }
        }
        else { // The start is contained, but the end is not : split the upper range
            
            //   |-------------------|         # Range to add
            //           |------------------|  # Other range whose boundary lies past the end of this one
            this.splitEndRange(newRange, containedEntry)
        }
        assert containedEntry.key < endPosition

        containedEntry.value.add(newRange)
    }
    
    @CompileStatic
    private void splitEndRange(IntRange newRange, Map.Entry<Integer,List<IntRange>> containedEntry) {
        
        final int endPosition = (int)newRange.to
        final int afterEndPosition = endPosition+1
        final List<IntRange> existingEndRanges = ranges[afterEndPosition]
        
        // Make a new range from the end of our range to the start of the next higher range
        // It needs to have the previous range only, not our new range
        List<IntRange> clonedList = new ArrayList(containedEntry.value.size())
        if(!existingEndRanges.is(null))
            clonedList.addAll(existingEndRanges)
            
        for(IntRange range in containedEntry.value) {
            if(endPosition < range.to )
                clonedList.add(range)
        }
        
        ranges[afterEndPosition] = clonedList
        checkRanges(afterEndPosition)
    }
    
    @CompileStatic
    private void splitStartRange(IntRange newRange, Map.Entry<Integer,List<IntRange>> lowerEntry) {
                
            final int startPosition = (int)newRange.from
            
            // Add all the overlapping regions from the adjacent lower region to
            // our new breakpoint
            List<IntRange> lowerSplitRegion =  new ArrayList(lowerEntry.value.size())
            for(IntRange lwrRange in lowerEntry.value) {
                if(lwrRange.to > startPosition)
                    lowerSplitRegion.add(lwrRange)
            }
                
            // Add our new range as covered by the split part
            lowerSplitRegion.add(newRange)
            ranges[startPosition] = lowerSplitRegion
            checkRanges(startPosition)
        
    }
    
    @CompileStatic
    private void addLowestRange(IntRange newRange) {
        
        final int startPosition = newRange.from
        final int endPosition = newRange.to
        final int positionAfterEnd = endPosition+1
        
        ranges.put(startPosition,[newRange])
        // checkRanges(startPosition)
        
        // If there are any ranges that start before the end of this new range,
        // add a breakpoint, and add this range to the intervening ranges
        int entryPos = startPosition
        Map.Entry<Integer, List<IntRange>> higherEntry = ranges.higherEntry(startPosition)
        Map.Entry<Integer, List<IntRange>> lastEntry = null
        while(higherEntry && higherEntry.key <= endPosition) {
            higherEntry.value.add(newRange)
            // checkRanges(higherEntry.key)
            lastEntry = higherEntry
            higherEntry = ranges.higherEntry(higherEntry.key)
        }
        
        // The last region needs to be split into two
        if(lastEntry) {
            if(!ranges.containsKey(positionAfterEnd)) {
              List<IntRange> newRanges = new ArrayList(lastEntry.value.size())
              for(IntRange r : lastEntry.value) {
                  if(positionAfterEnd>=r.from && positionAfterEnd<=r.to)
                      newRanges.add(r)
              }
              ranges[positionAfterEnd] = newRanges
              // checkRanges(endPosition+1)
            }
        }
        else {
            // Need to mark the ending of the range so that future adds will
            // look up the right point for adding new overlapping ranges
            if(!ranges.containsKey(positionAfterEnd)) {
                // println "Adding ${endPosition+1} to index ..."
                ranges[positionAfterEnd] = []
            }
        }
    }
    
    /**
     * Check that all the ranges added at a position span that position
     * 
     * @param position
     */
    @CompileStatic
    void checkRanges(int position) {
        return; // disabled!
        
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
    
    @CompileStatic
    List<IntRange> startingAt(int pos) {
        
        List<IntRange> result = this.ranges[pos]
        if(!result)
            return []
        
        // The given ranges could be starting at OR ending at the given position
        result.grep { IntRange r -> r.from == pos }
    }
    
    @CompileStatic
    List<IntRange> endingAt(int pos) {
        Map.Entry entry = this.ranges.lowerEntry(pos)
        if(!entry)
            return []
        List<IntRange> result = entry.value
        if(!result)
            return []
         result.grep { IntRange r -> r.to == pos }
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
    List<IntRange> getOverlaps(int start, int end) {
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
    
    final static Comparator<IntRange> INT_RANGE_COMPARATOR = new Comparator<IntRange>() {
            @Override
            @CompileStatic
            int compare(IntRange a, IntRange b) {
                if(a.is(b)) {
                    return 0
                }
                else {
                    int d = a.from - b.from; 
                    if(d == 0) {
                        d = a.to - b.to 
                        if(d == 0) {
                            d = System.identityHashCode(a).compareTo(System.identityHashCode(b))
                        }
                    }
                    return d
                }
            }
        } 
 
    
    /**
     * Return a list of the ranges in this index that overlap the given range.
     * If there are multiple ranges that overlap, or multiple distinct occurrences
     * of the same range that were added to the index, then these will be
     * returned separately.
     * 
     * @param start start of query interval 
     * @param end   end of query interval
     * @param returnFirst   if true, return only the first overlap found
     * 
     * @return List of ranges in the index that overlap the query
     */
    @CompileStatic
    List<Range> getOverlaps(final int start, final int end, final boolean returnFirst) {
        def resultSet = new TreeSet<IntRange>(INT_RANGE_COMPARATOR)
        def entry = ranges.lowerEntry(start+1)
        if(entry.is(null))
            entry = ranges.higherEntry(start)
            
        // Iterate as long as we are in the range
        while(!entry.is(null) && entry.key <= end) {
            for(IntRange r in entry.value) {
               if(r.from<=end && r.to>=start) {
                   resultSet.add(r)
                   if(returnFirst)
                       return (List<Range>)Arrays.asList(resultSet.toArray())
               }
            }
            entry = ranges.higherEntry(entry.key)
        }
        
        return (List<Range>)Arrays.asList(resultSet.toArray())
    }
    
    /**
     * Return a list of ranges that intersect the given start and end points.
     * <p>
     * <em>Note:</em>The start and end point are both treated as <b>inclusive</b>.
     * @param start start of range to find intersections for
     * @param end   end of range to find intersections for
     * @return  List of ranges that intersect the specified range.
     */
    @CompileStatic
    List<IntRange> intersect(int start, int end) {
       IntRange srcRange = start..end
       List<IntRange> result = getOverlaps(start,end)
       
       final int numResults = result.size()
       for(int i=0; i<numResults; ++i) {
           IntRange r = result[i]
           if(r instanceof GRange) 
               result[i] = (IntRange)((GRange)r).intersectRange(srcRange)
           else 
               result[i] = Math.max(r.from, start)..Math.min(r.to, end)            
       }
       return result
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
    @CompileStatic
    List<IntRange> subtractFrom(int start, int end) {

        List<IntRange> result = []
        List<IntRange> cutOut = intersect(start,end)       
        
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
            
        for(IntRange r in cutOut) {
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
    
    IntRange nearest(int pos) {
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
    
    /**
     * Returns the distance from the nearest range to the given position
     * 
     * @param pos
     * @return
     */
    int distanceTo(final int pos) {
        IntRange n = this.nearest(pos)
        if(n.containsWithinBounds(pos))
            return 0
        else
        if(n.from > pos) {
            return n.from - pos
        }
        else
        if(n.to < pos) {
            return pos - n.to
        }
        else
            assert false : "Position must be before, in or after the range"
    }
    
    @CompileStatic
    int distanceTo(final IntRange r) {
        int minToStart = distanceTo(r.from)
        int minToEnd = distanceTo(r.to)
        if(minToStart<0)
            return minToEnd
        else
        if(minToEnd<0)
            return minToStart
        else
            return Math.min(minToEnd,minToStart)
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
    
    void debugdump() {
        this.ranges.each { int pos, List intervalRanges ->
            println "$pos\t:" + intervalRanges.collect { it.from + '-' + it.to }.join(",")
        }
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
    @CompileStatic
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
            
            @CompileStatic
            boolean hasNext() {
                if(nextRange)
                    return true
                    
                // Note: pos becomes null when iteration through index exhausted
                // nextRange is set null at each call of next()
                while(!pos.is(null) && (nextRange.is(null) || nextRange.from < pos)) {
                  nextRange = findNext()
                }
                return !nextRange.is(null)
            }
            
            @CompileStatic
            IntRange next() {
                
                if(!hasNext()) // Note: populates nextRange
                    throw new IndexOutOfBoundsException("No more ranges in this index")
                
                IntRange result = nextRange
                nextRange = null
                return result
            }
            
            @CompileStatic
            IntRange findNext() {
                if(pos == -1) {
                    if(ranges.isEmpty())
                        pos = null
                    else
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
     * representing all the ranges covered by any range in this RangeIndex.
     * 
     * @return
     */
    @CompileStatic
    RangeIndex reduce(Closure reducer = null) {
        RangeIndex reduced = new RangeIndex()

        // We take advantage of the fact that we iterate the ranges in order of genomic start position
        IntRange currentRange = null
        for(IntRange r in this) {
            if(null == currentRange) {
                currentRange = r
                continue
            }
            
            assert r.from >= currentRange.from
            
            if(GRange.overlaps(r,currentRange)) { 
                if(currentRange instanceof GRange) {
                    Object newExtra = currentRange.getExtra()
                    if(r instanceof GRange) {
                        if(newExtra == null)
                            newExtra = ((GRange)r).getExtra()
                        else
                        if(reducer != null) { // Both are GRanges and both have extra values - they need to be merged
                            newExtra = reducer(currentRange, r)
                        }
                    }
                    
                    currentRange = new GRange(Math.min((int)r.from, currentRange.from),Math.max((int)r.to, currentRange.to), null)
                    
                    // Regions as extras violate the assumption that a Region extra reflects 
                    // the same region as a range. So we create a new region, but we 
                    // preserve the properties
                    if(newExtra instanceof Region) {
                        Map props = ((Expando)newExtra).getProperties()
                        newExtra = new Region(props,((Region)newExtra).chr, currentRange)
                    }
                    currentRange.extra = newExtra
                }
                else {
                    Object newExtra = r instanceof GRange ? ((GRange)r).extra : null
                    currentRange = new GRange(Math.min((int)r.from, currentRange.from),Math.max((int)r.to, currentRange.to),newExtra)
                }
            }
            else {
                reduced.add(currentRange)
                currentRange = r
            }
        }
        
        if(null != currentRange) {
            reduced.add(currentRange)
        }
            
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
            
            Map.Entry<Integer,List<IntRange>> prev = null
            
            IntRange prevRange = null
            
            boolean hasNext() {
//                println "hasNext = " + i.hasNext() + " prev = " + prev
                i.hasNext() || (prev.is(null) && !prev.value.empty)
            }
            
            GRange next() {
                
                Map.Entry<Integer,List<IntRange>> curr = null
                if(prev.is(null)) {
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
                    if(!prevRange.is(null)) {
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
    
    /**
     * Return a list of ranges representing "coverage" blocks within this
     * index. That is, for each position where the number of overlapping
     * ranges changes, a separate range is returned with a count of the overlaps
     * as the 'extra' field.
     *  
     * @return  list of ranges representing coverage blocks with number of overlaps
     *          as the 'extra' field.
     */
    List<GRange> coverage() {
        
        List<GRange> result = []
        
        Iterator<Map.Entry<Integer,List<IntRange>>> iter = ranges.iterator()
        GRange prevRange = null
        List<IntRange> prevRanges
        for(Map.Entry<Integer,List> rangesEntry in iter) {
            
            List<IntRange> ranges = rangesEntry.value
            GRange covRange = null
            if(prevRange != null) {
                covRange = new GRange(prevRange.to,  rangesEntry.key, prevRanges.size())
                result << covRange
            }
            else
                covRange = new GRange(0, rangesEntry.key, null)
                
            prevRanges = ranges
            prevRange = covRange    
        }
        return result
    }
}