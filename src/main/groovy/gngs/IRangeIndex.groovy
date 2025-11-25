package gngs

import groovy.transform.CompileStatic

/**
 * Common interface for range indexing implementations.
 * Defines methods for adding, removing, and querying ranges in an index.
 */
/**
 * Common interface for range indexing implementations.
 * Defines methods for adding, removing, and querying ranges in an index.
 */
interface IRangeIndex extends Iterable<IntRange> {
    
    /**
     * Add the given range to this index. The end coordinate is treated as exclusive
     * to match BED convention and RangeIndex behavior.
     */
    void add(IntRange range)
    
    /**
     * Add range with explicit start/end coordinates
     * If the 'extra' argument is specified then a GRange will be added with the extra
     * as its extra data argument. Otherwise, a plain IntRange will be added.
     */
    void add(int start, int end, Object extra)
    
    /**
     * Remove an existing range from the index. The existing range
     * MUST already be in the index, and cannot be a user-supplied range.
     */
    void remove(IntRange r)
    
    /**
     * Get the number of ranges in the index
     */
    int getNumRanges()
    
    /**
     * Get ranges that start at the given position
     */
    List<IntRange> startingAt(int pos)
    
    /**
     * Get ranges that end at the given position (exclusive)
     */
    List<IntRange> endingAt(int pos)
    
    /**
     * Find all ranges that overlap the given position
     */
    List<IntRange> getOverlaps(int position)
    
    /**
     * Return a list of ranges stored in the index that overlap the specified
     * range. Both ends of the range are *inclusive*.
     * 
     * @param start first position to look for overlaps
     * @param end   last position to look for overlaps (inclusive)
     * @return  List of ranges overlapping the specified start -> end range
     */
    List<IntRange> getOverlaps(int start, int end)
    
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
    List<IntRange> getOverlaps(int start, int end, boolean returnFirst)
    
    /**
     * Return true if the given range overlaps at least one range in this
     * RangeIndex.
     * 
     * @param start start of range to test (inclusive)
     * @param end   end point of range to test (inclusive)
     * @return  true iff at least one range in this index overlaps the given start and end point
     */
    boolean overlaps(int start, int end)
    
    /**
     * Return a list of ranges that intersect the given start and end points.
     * <p>
     * <em>Note:</em>The start and end point are both treated as <b>inclusive</b>.
     * @param start start of range to find intersections for
     * @param end   end of range to find intersections for
     * @return  List of ranges that intersect the specified range.
     */
    List<IntRange> intersect(int start, int end)
    
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
    List<IntRange> subtractFrom(int start, int end)
    
    /**
     * Find the next range after the given position
     */
    Range nextRange(int pos)
    
    /**
     * Find the previous range before the given position
     */
    Range previousRange(int pos)
    
    /**
     * Find the nearest range to the given position
     */
    IntRange nearest(int pos)
    
    /**
     * Returns the distance from the nearest range to the given position
     */
    int distanceTo(int pos)
    
    /**
     * Returns the distance from this range index to the given range
     */
    int distanceTo(IntRange r)
    
    /**
     * Get the first range in the index
     */
    Range first()
    
    /**
     * Get the last range in the index
     */
    Range last()
    
    /**
     * Return an iterator that will return each range in the index
     * in genomic order by starting position, starting from the given 
     * position.
     */
    Iterator<IntRange> iteratorAt(int startingPos)
    
    /**
     * Return an iterator over the entire tree that returns intervals in reverse order.
     */
    Iterator<IntRange> reverseIterator()
    
    /**
     * Return an iterator that will return each range that starts
     * at or before the given position, proceeding backwards
     * through the genome, ordered by start position.
     */
    Iterator<IntRange> reverseIteratorAt(int startingPos)
    
    /**
     * Merge all overlapping ranges together to make simplified regions
     * representing all the ranges covered by any range in this RangeIndex.
     * 
     * If a closure is provided, it is called whenever two regions need to be merged,
     * providing both regions as arguments. The return value is set as the <code>extra</code>
     * attribute on the resulting region.
     * 
     * @return  new IRangeIndex with all overlapping regions reduced to single flattened regions
     */
    IRangeIndex reduce(Closure reducer)
    
    /**
     * Return a list of ranges representing "coverage" blocks within this
     * index. That is, for each position where the number of overlapping
     * ranges changes, a separate range is returned with a count of the overlaps
     * as the 'extra' field.
     *  
     * @return  list of ranges representing coverage blocks with number of overlaps
     *          as the 'extra' field.
     */
    List<GRange> coverage()
}
