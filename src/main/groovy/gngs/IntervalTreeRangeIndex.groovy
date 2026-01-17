package gngs

import groovy.transform.CompileStatic
import htsjdk.samtools.util.IntervalTree

/**
 * An implementation of RangeIndex backed by IntervalTree for more efficient
 * overlap queries. This implementation uses inclusive end coordinates internally
 * but maintains the same BED-style (exclusive end) interface as RangeIndex.
 */
@CompileStatic
class IntervalTreeRangeIndex implements IRangeIndex {
    
    private IntervalTree<List<IntRange>> tree = new IntervalTree<>()
    
    private int numRanges = 0
    
    /**
     * Add the given range to this index. The end coordinate is treated as exclusive
     * to match BED convention and RangeIndex behavior.
     */
    void add(IntRange range) {
        // Convert to inclusive end for IntervalTree
        // For single-base intervals, use same position for start and end
        int endPos = range.from == range.to ? range.from : range.to-1
        def existing = tree.find(range.from, endPos)?.getValue()
        if (existing) {
            existing.add(range)
        } else {
            tree.put(range.from, endPos, [range])
        }
        ++numRanges
    }
    
    /**
     * Add range with explicit start/end coordinates
     */
    void add(int start, int end, Object extra = null) {
        IntRange range = extra != null ? new GRange(start, end-1, extra) : new IntRange(start, end-1)
        add(range)
    }
    
    /**
     * Get the number of ranges in the index
     */
    int getNumRanges() {
        return numRanges
    }
    
    /**
     * Find all ranges that overlap the given position
     */
    List<IntRange> getOverlaps(int position) {
        return tree.overlappers(position, position)
            .collect { it.getValue() }
            .collectMany { it } as List<IntRange>
    }
    
    /**
     * Find all ranges that overlap the given range. End coordinate is exclusive.
     */
    List<IntRange> getOverlaps(int start, int end) {
        return getOverlaps(start, end, false)
    }
    
    /**
     * Find overlaps with option to return after first match
     */
    List<IntRange> getOverlaps(int start, int end, boolean returnFirst) {
        // For single-base queries, use same position for start and end
        int endPos = start == end ? start : end
        def overlaps = tree.overlappers(start-1, endPos+1)
        
        def intersects = overlaps.collectMany { 
                it.getValue() 
            }
            .findAll { range -> 
                // A range overlaps if:
                // - it starts at or before the query end AND
                // - it ends at or after the query start
                // Note: range.to is exclusive, so we need to subtract 1 for comparison
                range.from <= end && start <= range.to
            } as List<IntRange>

        return returnFirst ? intersects.take(1) : intersects
    }
    
    /**
     * Test if any range overlaps the given coordinates
     */
    boolean overlaps(int start, int end) {
        !getOverlaps(start, end, true).isEmpty()
    }
    
    /**
     * Remove a range from the index
     */
    void remove(IntRange r) {
        // Convert to inclusive end for IntervalTree
        tree.remove(r.from, r.to-1)
        --numRanges
    }
    
    /**
     * Iterate over all ranges in order
     */
    @Override
    Iterator<IntRange> iterator() {
        return tree.iterator()
            .collect { it.getValue() }
            .collectMany { it }
            .iterator() as Iterator<IntRange>
    }
    
    /**
     * Get ranges that start at the given position
     */
    List<IntRange> startingAt(int pos) {
        return tree.find(pos, pos)?.getValue() ?: []
    }
    
    /**
     * Get ranges that end at the given position (exclusive)
     */
    List<IntRange> endingAt(int pos) {
        // Need to check pos-1 since we store inclusive ends
        return tree.find(pos-1, pos-1)?.getValue() ?: []
    }
    
    /**
     * Support for Groovy 'in' operator by implementing isCase
     * Tests if the given position falls within any range in this index
     */
    boolean isCase(int position) {
        return tree.overlappers(position, position).hasNext()
    }
    
    /**
     * Merge all overlapping ranges together to make simplified regions
     * representing all the ranges covered by any range in this RangeIndex.
     * 
     * If a closure is provided, it is called whenever two regions need to be merged,
     * providing both regions as arguments. The return value is set as the <code>extra</code>
     * attribute on the resulting region.
     * 
     * @return new IRangeIndex with all overlapping regions reduced to single flattened regions
     */
    @CompileStatic
    IRangeIndex reduce(Closure reducer = null) {
        IRangeIndex reduced = new IntervalTreeRangeIndex()
        
        // No ranges to reduce
        if(numRanges == 0) {
            return reduced
        }
        
        // Process ranges in start position order
        Iterator<IntervalTree.Node<List<IntRange>>> iter = tree.iterator()
        IntRange currentRange = null
        
        while(iter.hasNext()) {
            List<IntRange> ranges = iter.next().getValue()
            for(IntRange r : ranges) {
                if(currentRange == null) {
                    currentRange = r
                    continue
                }
                
                // If current range overlaps with this range
                if(r.from <= currentRange.to) {
                    // Merge the ranges
                    if(currentRange instanceof GRange) {
                        Object newExtra = ((GRange)currentRange).extra
                        if(r instanceof GRange) {
                            if(reducer != null) {
                                newExtra = reducer(currentRange, r)
                            }
                        }
                        currentRange = new GRange(
                            Math.min(currentRange.from, r.from),
                            Math.max(currentRange.to, r.to),
                            newExtra
                        )
                    }
                    else {
                        currentRange = new IntRange(
                            Math.min(currentRange.from, r.from),
                            Math.max(currentRange.to, r.to)
                        )
                    }
                }
                else {
                    // No overlap - add current range and start new one
                    reduced.add(currentRange)
                    currentRange = r
                }
            }
        }
        
        // Add final range if exists
        if(currentRange != null) {
            reduced.add(currentRange)
        }
        
        return reduced
    }
    
    /**
     * Not yet implemented
     */
    List<GRange> coverage() {
        throw new UnsupportedOperationException("coverage() not yet implemented") 
    }
    
    /**
     * Return a list of ranges that intersect the given start and end points.
     * <p>
     * <em>Note:</em>The start and end point are both treated as <b>inclusive</b>.
     */
    List<IntRange> intersect(int start, int end) {
        def overlaps = getOverlaps(start, end)
        return overlaps.collect { range ->
            // For each overlapping range, create a new range representing 
            // the intersection with the query range
            int intersectStart = Math.max(start, range.from)
            int intersectEnd = Math.min(end, range.to)
            new IntRange(intersectStart, intersectEnd)
        }
    }
    
    /**
     * Subtract all the ranges in this range index from the given range and return
     * the resulting list of ranges.
     * <p>
     * <em>Note:</em>both start and end are considered <b>inclusive</b>.
     * 
     * @param start start of range to subtract from (inclusive)
     * @param end   end of range to subtract from (inclusive)
     * @return  List of ranges left after the ranges in this RangeIndex are removed from 
     *          the specified region.
     */
    List<IntRange> subtractFrom(int start, int end) {
        // Get all overlapping ranges, sorted by start position
        List<IntRange> overlaps = getOverlaps(start, end).sort { it.from }
        
        // If no overlaps, return the entire range
        if (overlaps.isEmpty()) {
            return [new IntRange(start, end)]
        }
        
        List<IntRange> result = []
        int currentPos = start
        
        for (IntRange overlap : overlaps) {
            // If there's a gap before this overlap, add it to result
            if (currentPos < overlap.from) {
                result.add(new IntRange(currentPos, overlap.from - 1))
            }
            
            // Move current position past this overlap
            // overlap.to is exclusive in our representation, so we use it directly
            currentPos = Math.max(currentPos, overlap.to + 1)
            
            // If we've gone past the end, we're done
            if (currentPos > end) {
                break
            }
        }
        
        // If there's remaining range after all overlaps, add it
        if (currentPos <= end) {
            result.add(new IntRange(currentPos, end))
        }
        
        return result
    }
    
    /**
     * Find the next range after the given position
     */
    Range nextRange(int pos) {
        // Use IntervalTree's min() method to efficiently find the first range
        // that starts after the given position
        IntervalTree.Node<List<IntRange>> node = tree.min(pos + 1, pos + 1)
        
        if (node == null) {
            return null
        }
        
        // The node may contain multiple ranges at the same position
        // Return the first one (they all have the same start position)
        List<IntRange> ranges = node.getValue()
        return ranges.isEmpty() ? null : ranges[0]
    }
    
    /**
     * Find the previous range before the given position.
     * Returns the range whose end point is before the given position.
     */
    Range previousRange(int pos) {
        // Find the first range that ends before pos by searching backwards
        // We look for ranges that could potentially end before pos
        IntRange candidate = null
        
        // Iterate backwards to find a range that ends before pos
        Iterator<IntervalTree.Node<List<IntRange>>> iter = tree.reverseIterator(pos - 1, pos - 1)
        while (iter.hasNext()) {
            List<IntRange> ranges = iter.next().getValue()
            for (IntRange range : ranges) {
                // Check if this range ends before pos (remember: range.to is exclusive)
                if (range.to < pos) {
                    candidate = range
                    break
                }
            }
            if (candidate != null) {
                break
            }
        }
        
        // If no candidate found, there's no range ending before pos
        if (candidate == null) {
            return null
        }
        
        // Now check if there are any ranges that overlap the candidate's endpoint
        // and extend closer to pos. If so, the one with the highest endpoint is the answer.
        List<IntRange> overlappers = getOverlaps(candidate.to, pos - 1)
        
        // If no overlappers, candidate is the answer
        if (overlappers.isEmpty()) {
            return candidate
        }
        
        // Find the overlapper with the highest end point that's still < pos
        IntRange best = candidate
        for (IntRange overlapper : overlappers) {
            if (overlapper.to < pos && overlapper.to > best.to) {
                best = overlapper
            }
        }
        
        return best
    }
    
    /**
     * Find the nearest range to the given position
     */
    IntRange nearest(int pos) {
        throw new UnsupportedOperationException("nearest() not yet implemented")
    }
    
    /**
     * Returns the distance from the nearest range to the given position
     */
    int distanceTo(int pos) {
        throw new UnsupportedOperationException("distanceTo(int) not yet implemented")
    }
    
    /**
     * Returns the distance from this range index to the given range
     */
    int distanceTo(IntRange r) {
        throw new UnsupportedOperationException("distanceTo(IntRange) not yet implemented")
    }
    
    /**
     * Get the first range in the index
     */
    Range first() {
        throw new UnsupportedOperationException("first() not yet implemented")
    }
    
    /**
     * Get the last range in the index
     */
    Range last() {
        throw new UnsupportedOperationException("last() not yet implemented")
    }
    
    /**
     * Return an iterator that will return each range in the index
     * in genomic order by starting position, starting from the given 
     * position.
     */
    Iterator<IntRange> iteratorAt(int startingPos) {
        throw new UnsupportedOperationException("iteratorAt() not yet implemented")
    }
    
    /**
     * Return an iterator over the entire tree that returns intervals in reverse order.
     */
    Iterator<IntRange> reverseIterator() {
        throw new UnsupportedOperationException("reverseIterator() not yet implemented")
    }
    
    /**
     * Return an iterator that will return each range that starts
     * at or before the given position, proceeding backwards
     * through the genome, ordered by start position.
     */
    Iterator<IntRange> reverseIteratorAt(int startingPos) {
        throw new UnsupportedOperationException("reverseIteratorAt() not yet implemented")
    }
}
