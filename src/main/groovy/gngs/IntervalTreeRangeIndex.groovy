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
        def overlaps = tree.overlappers(start, endPos)
            .collect { it.getValue() }
            .collectMany { it }
            .findAll { range -> 
                // A range overlaps if:
                // - it starts at or before the query end AND
                // - it ends at or after the query start
                range.from <= end && start <= range.to
            } as List<IntRange>
        return returnFirst ? overlaps.take(1) : overlaps
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
     * Not yet implemented
     */
    RangeIndex reduce(Closure reducer = null) {
        throw new UnsupportedOperationException("reduce() not yet implemented")
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
        throw new UnsupportedOperationException("intersect() not yet implemented")
    }
    
    /**
     * Subtract all the ranges in this range index from the given range and return
     * the resulting list of ranges.
     */
    List<IntRange> subtractFrom(int start, int end) {
        throw new UnsupportedOperationException("subtractFrom() not yet implemented")
    }
    
    /**
     * Find the next range after the given position
     */
    Range nextRange(int pos) {
        throw new UnsupportedOperationException("nextRange() not yet implemented")
    }
    
    /**
     * Find the previous range before the given position
     */
    Range previousRange(int pos) {
        throw new UnsupportedOperationException("previousRange() not yet implemented")
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
