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
 * A set of genomic regions, each indexed by chromosome, start and end position.
 * <p>
 * The Regions class provides operations for working on Ranges but adds awareness of chromosomes 
 * (or contigs) to {@link RangeIndex}. It offers most of the same functions as exist for a 
 * {@link RangeIndex} but transparently handles multiple chromosomes.
 * <p>
 * The internal representation for each genomic range stored in a {@link Regions} object
 * is a Groovy {@link IntRange}. For efficiency, many of the methods return these
 * underlying IntRange objects directly, whenever a method is clearly appliciable only
 * to a single chromosome. When functionality spans chromosomes, the {@link IntRange} 
 * objects will be converted to {@link Region} objects which carries an extra field for 
 * the chromosome. This adds a small amount of overhead, but enables the class to 
 * support various useful functions such as implementing Groovy's enhanced support for
 * {@link Iterable} objects.
 * <p>
 * Often Regions will be loaded from a BED file. In that case, see the {@link BED} class
 * which extends Regions, but provides extensive support for parsing, loading and 
 * filtering BED files. To create a purely in-memory set of regions, use the constructor:
 * <pre>
 * Regions regions = new Regions()
 * regions.addRegion("chr1",100,200)
 * regions.addRegion("chr1",150,300)
 * </pre>
 * When added this way, Regions follows the BED file convention that the start of the range
 * is considered <i>inclusive</i>, but the end value is considered <i>exclusive</i> of 
 * the range added. Once loaded, various operations can be performed on the ranges. Since
 * the Iterable interface is implemented, you can use any of the Groovy special operations
 * that work on iterables:
 * <pre>
 * int bases = regions.grep { it.chr == "chrX" || it.chr == "chrY" }*.size().sum()
 * println "There were $bases bases from sex chromosomes"
 * </pre>
 * The Regions class offers many operations for convenient querying of and logical 
 * operations on intervals:
 * <b>Finding overlaps:</b>
 * <pre>
 * Region r = new Region("chr1",120,130)
 * regions.getOverlaps(r).size()==1
 * </pre>
 * <p>Many operations are possible using built in Groovy iterator methods. for example, to find 
 * the indexes of overlapping regions:</p>
 * <pre>
 * assert regions.findIndexValues { r.overlaps(it) } == [ 0 ]
 * </pre>
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class Regions implements Iterable<Region> {
        
    /**
     * Index for looking up overlaps
     */
    Map<String, RangeIndex> index = [:]
    
    /**
     * A list of ranges in the order they were loaded
     */
    Map<String,List<IntRange>> allRanges = [:]
    
    /**
     * Create new empty set of regions 
     */
    Regions(Map attributes=[:]) {
    }
    
    /**
     * Create a Regions from a list of Ranges
     */
    Regions(Map attributes=[:], String chr, Iterable<Range> ranges) {
        for(Range r in ranges) {
            addRegion(chr, r.from, r.to)
        }
    }
    
    /**
     * Create a Regions from a list of regions
     */
    Regions(Map attributes=[:], Iterable<IRegion> regions) {
        for(IRegion r in regions) {
            
            if(r instanceof Region) {
                this.addRegion((Region)r)
                continue
            }
            
            Range range = r.range
            if(r.range instanceof GRange)
                addRegion(r.chr, range.from, range.to, range.extra)
            else
                addRegion(r.chr, range.from, range.to)
        }
    }
    
    /**
     * Return a new Regions object that has <code>bp</code> bases
     * added to the beginning and end of each interval. If this
     * causes overlaps then these are left in the resulting 
     * Objects.
     * 
     * @param bp
     * @return  a new Regions with bp wider intervals
     */
    @CompileStatic
    Regions widen(int bp) {
        Regions regions = new Regions()
        for(Region r in this) {
            Region widened = r.widen(bp)
            if(!r.properties.is(null) && !r.properties.isEmpty()) {
                widened.setProps(r.properties)
            }
            regions.addRegion(widened)
        }
        return regions
    }
    
    @CompileStatic
    long size() {
        long sizeCount = 0
        this.eachRange { String chr, int start, int end ->
            sizeCount += (end - start + 1)
        }
        return sizeCount
    }
    
    @CompileStatic
    void eachRange(Closure c) {
       eachRange(unique:false, c)
    }
    
    /**
     * Iterate over each region and invoke the specified
     * closure with bed file information for the line.
     * <p>
     * Options can be provided as a Map argument, including:
     *
     * <li>unique (whether identical ranges should be included multiple times or only once)
     */
//    @CompileStatic
    void eachRange(Map options, Closure c) {
        ProgressCounter progress = new ProgressCounter(lineInterval:10, timeInterval:10000)
        boolean includeInfo = c.getMaximumNumberOfParameters()==4;
        boolean useRange = c.getMaximumNumberOfParameters()==2;
        boolean useRegion = c.getMaximumNumberOfParameters()==1;
        boolean unique = options.unique?true:false
        this.allRanges.each { String chr, List<groovy.lang.IntRange> ranges ->
            HashSet processed = new HashSet()
            for(groovy.lang.IntRange r in ranges) {
                
              if(unique) {
                  String key = r.from + ":" + r.to
                  if(processed.contains(key)) {
                      progress.count()
                      continue
                  }
                  processed.add(key)
              }
  
              if(includeInfo) {
                  c.call(chr,r.from,r.to,r instanceof GRange ? ((GRange)r).getExtra() : null)
              }
              else
              if(useRange) {
                  c.call(chr,r)
              }
              else
              if(useRegion) {
                  c.call(new Region(chr, r))
              }
              else {
                  c.call(chr,r.from,r.to)
              }
                  
              progress.count()
            }
        }
    }
    
    @CompileStatic
    List<IntRange> intersect(IRegion region) {
        intersect(region.chr, (int)region.range.from, (int)region.range.to)
    }
    
    @CompileStatic
    List<IntRange> intersect(String chr, int start, int end) {
        RangeIndex chrIndex = this.index[chr]
        if(chrIndex == null)
            return []
        List<IntRange> result = chrIndex.intersect(start,end)
        return result
    }
    
    @CompileStatic
    Regions intersectRegions(Regions other) {
        this.intersect(other)
    }
    
    @CompileStatic
    Regions intersectRegion(Region other) {
        Regions result = new Regions()
        for(r in intersect(other)) {
            result.addRegion(new Region(other.chr, r.from, r.to))
        }
        return result
    }
    
    @CompileStatic
    Regions intersect(BED other) {
        this.intersectImpl(other)
    }
    
    @CompileStatic
    Regions intersect(Regions other) { 
        this.intersectImpl(other)
    }
    
    /**
     * Returns a set of regions representing each region in this Regions
     * intersected with the Regions in the other regions.
     * <p>
     * Note: nothing is done to deal specially with overlapping ranges.
     * Thus if there are overlapping ranges in this Regions then you will
     * find overlapping ranges in the results wherever those overlaps intersect
     * with the other Regions. The same applies with respect to the other Regions.
     * For a "flat" intersection of the two Regions, you should use {@link #reduce()} to
     * flatten the source and target first before calling this method.
     */
    @CompileStatic
    Regions intersectImpl(Regions other) {
        Regions result = new Regions()
        this.index.each { chr, chrIndex ->
            RangeIndex otherChrIndex = other.index[chr]
            if(otherChrIndex.is(null))
                return 
            for(IntRange r in chrIndex) {
                List<IntRange> otherRanges = otherChrIndex.intersect(r.from, r.to)
                for(IntRange xr : otherRanges) { result.addRegion(chr, (int)xr.from, ((int)xr.to)+1) }
            }
        }
        return result
    }
    
    /**
     * Find the ranges that have at least 1bp overlap with the given region
     * <p>
     * Note that this method returns a list of {@link groovy.lang.IntRange} objects.
     * Where these belong to Region objects they will be instances of GRange objects
     * with the Region set as the {@link gngs.Region#extra} field.
     * <p>
     * If you want to get the {@link gngs.Region} objects back directly, use
     * {@link gngs.Region#getOverlapRegions()}.
     * 
     * @see         gngs.Region#getOverlapRegions()
     * @param       the region to locate overlaps for0
     */
    @CompileStatic
    List<IntRange> getOverlaps(IRegion r) {
        getOverlaps(r.chr, r.range.from, r.range.to)
    }
    
    /**
     * Returns the overlaps with with the given region as region objects.
     * <p>
     * Note: all the internal ranges must be stored as full Region objects,
     * or this operation will throw an exception.
     */
    @CompileStatic
    List<Region> getOverlapRegions(IRegion r) {
        getOverlaps(r.chr, r.range.from, r.range.to).collect { IntRange overlap ->
            (Region)((GRange)overlap).extra
        }
    }
    
    /**
     * Return true if the given region overlaps any range in this Regions
     * 
     * @param r Region to test for overlaps
     * @return  true iff the region overlaps at least one region in this Regions
     */
    @CompileStatic
    boolean overlaps(String chr, int from, int to) {
        RangeIndex chrIndex = this.index[chr]
        if(chrIndex == null)
            return false
        return chrIndex.overlaps(from,to)        
    }
    
    
    /**
     * Return true if the given region overlaps any range in this Regions
     * 
     * @param r Region to test for overlaps
     * @return  true iff the region overlaps at least one region in this Regions
     */
    @CompileStatic
    boolean overlaps(IRegion r) {
        return this.overlaps(r.chr, r.range.from,r.range.to)        
    }
    
    /**
     * Returns true if at least one region in r overlaps at least one
     * region in this Regions object
     * 
     * @param r
     * @return
     */
    @CompileStatic
    boolean overlaps(Regions other) {
        for(Region r in other ) {
            if(this.overlaps(r))
                return true
        }
        return false
    }
    
    /**
     * Returns true if at least one region in r overlaps at least one
     * region in this Regions object
     * 
     * @param r
     * @return
     */
    @CompileStatic
    boolean overlaps(Iterable<IRegion> other) {
        for(IRegion r : other) {
            if(this.overlaps(r))
                return true
        }
        return false
    }
     
    /**
     * Return a list of ranges that overlap the specified range. 
     * <em>Note</em>: both ends of the range are *inclusive*.
     * The Range objects returned all belong to the reference sequence <em>chr</em>.
     * Note, if the internal ranges are GRanges objects (as they will be often)
     * then you can get full Region objects out from the 'extra' property:
     * <pre>
     * List<Region> overlaps = regions.getOverlaps('chrX',10000,20000)*.extra
     *</pre>
     * @param start first position to look for overlaps
     * @param end   last position to look for overlaps (inclusive)
     * @return  List of ranges overlapping the specified start -> end range
     */
    @CompileStatic
    List<IntRange> getOverlaps(String chr, int start, int end) {
        RangeIndex chrIndex = this.index[chr]
        if(chrIndex == null)
            return []
        return chrIndex.getOverlaps(start,end)
    }
    
    @CompileStatic
    List<Region> subtractFrom(Region region) {
        subtractFrom(region.chr, region.from, region.to+1).collect {
            new Region(region.chr, it.from..it.to)
        }
    }
    
    /**
     * Remove all of the regions belonging to this {@link Regions} object
     * from the interval specified, and return the list of Ranges that results.
     * <p>
     * Note: the end attribute is considered <i>exclusive</i> of the range to be
     * subtracted from.
     * 
     * @param chr   chromsome of range
     * @param start start of interval to subtract from (inclusive)
     * @param end   end of interval to subtract from (exclusive)
     * @return
     */
    @CompileStatic
    List<IntRange> subtractFrom(String chr, int start, int end) {
        RangeIndex chrIndex = this.index[chr]
        if(chrIndex == null)
            return [start..end-1]
        
        return chrIndex.subtractFrom(start,end+1)
    }
    
    /**
     * Call closure c for each range that overlaps the given position
     */
    @CompileStatic
    void eachOverlap(String chr, int pos, Closure c) {
        RangeIndex chrIndex = this.index[chr]
        if(chrIndex == null)
            return
        for(groovy.lang.Range r in chrIndex.getOverlaps(pos)) {
            c(chr, r.from, r.to)
        }
    }
    
    @CompileStatic
    List<Region> regionsStartingAt(String chr, int pos) {
        return (List<Region>)((List<GRange>)startingAt(chr, pos))*.extra
    }
    
    @CompileStatic
    List<IntRange> startingAt(String chr, int pos) {
        RangeIndex chrIndex = this.index[chr]
        if(chrIndex == null)
            return []
        return chrIndex.startingAt(pos)
    }
    
    /**
     * Return a list of ranges that end exactly at the specified position.
     *
     * @param chr
     * @param pos
     * @return
     */
    List<Region> regionsEndingAt(String chr, int pos) {
        return endingAt(chr, pos)*.extra
    } 
    
    /**
     * Return a list of ranges that end exactly at the specified position.
     * <p>
     * NOTE: the position is considered <i>inclusive</i> to the range. This
     * is different to BED file notation.
     *
     * @param chr
     * @param pos
     * @return
     */
    List<Range> endingAt(String chr, int pos) {
        RangeIndex chrIndex = this.index[chr]
        if(chrIndex == null)
            return []
        
        return chrIndex.endingAt(pos)
    }
    
    /**
     * Returns the prior range that has its end closest
     * to the given position.
     */
    groovy.lang.IntRange previousRange(String chr,int pos) {
        RangeIndex chrIndex = this.index[chr]
        if(chrIndex == null)
            return null
        return chrIndex.previousRange(pos)
    }
    
    /**
     * Returns the next range that has its beginning
     * closest to the given position.
     */
    @CompileStatic
    groovy.lang.IntRange nextRange(String chr, int pos) {
        RangeIndex chrIndex = this.index[chr]
        if(chrIndex == null)
            return null
        return (IntRange)chrIndex.nextRange(pos)
    }
    
    /**
     * Returns the region that is count regions forward of the given position
     */
    groovy.lang.IntRange forward(String chr, int pos, int count = 1) {
        RangeIndex chrIndex = this.index[chr]
        if(chrIndex == null)
            return null
            
        groovy.lang.IntRange result
        Iterator<IntRange> iter = chrIndex.iteratorAt(pos)
        while(iter.hasNext() && count > 0) {
            result = iter.next()
            --count
        }
        return result
    }
    
    /**
     * Returns the region that is count regions backwards from the given position
     */
    groovy.lang.IntRange backward(String chr, int pos, int count = 1) {
        RangeIndex chrIndex = this.index[chr]
        if(chrIndex == null)
            return null
            
        groovy.lang.IntRange result
        Iterator<IntRange> iter = chrIndex.reverseIteratorAt(pos)
        while(iter.hasNext() && count > 0) {
            result = iter.next()
            
            // note: when pos is the start of a region, the first region
            // returned from the reverse iterator is that region - but this does
            // not represent going backward at all, so only count if the 
            // region returned starts before pos
            if(result.from < pos)
                --count
        }
        return result
    }

    
    /**
     * Return a window of n regions upstream and downstream of the 
     * given region
     */
    List<Region> window(Region r, int n) {
        
        List<Region> result = []
        
        RangeIndex index = this.index[r.chr]
        Iterator downstream = index.reverseIteratorAt(r.from)
        int downstreamCount = 0
        while(downstream.hasNext() && downstreamCount <= n) {
            Region downR = new Region(r.chr, downstream.next())
            if(downR.to != r.to || downR.from != r.from) {
                ++downstreamCount
                result.add(0, downR)
            }
        }
        
        result << r
        
        int upstreamCount = 0
        Iterator upstream = index.iteratorAt(r.to)
        while(upstream.hasNext() && upstreamCount < n) {
            Region upR = new Region(r.chr, upstream.next())
            if(upR.to != r.to || upR.from != r.from) {
                result << upR
                ++upstreamCount
            }
        }
        return result
    }
    
    @CompileStatic
    void remove(String chr, Range r) {
        RangeIndex chrIndex = this.index[chr]
        if(!chrIndex)
            throw new IllegalArgumentException("Cannot remove region from chrosomome $chr : it does not have any regions")
        chrIndex.remove(r)
        
        if(this.allRanges[chr]) // Could be really slow?
            this.allRanges[chr].removeAll { Range r2 -> r.from == r2.from && r.to == r2.to }
    }
    
    /**
     * Returns the range that is "closest" to the given
     * position.
     *
     *   1. if no range overlaps, the closer of a) the end
     *      of the nearest prior range vs b) the start of the
     *      nearest following range is returned.
     *
     *   2. if one or more ranges overlap, returns the overlapping
     *      range with the start or end closest to the specified position
     *
     * @param chr
     * @param pos
     * @return
     */
    groovy.lang.IntRange nearest(String chr, int pos) {
        this.index[chr]?.nearest(pos)
    }
    
    /**
     * @return the distance to the nearest region in this Regions object, or -1 
     *         if there is no region in the same contig, if a region onverlaps, 
     *         returns 0
     */
    @CompileStatic
    int distanceTo(String chr, int pos) { 
        return this.index.containsKey(chr) ? this.index[chr].distanceTo(pos) : (int)-1
    }
    
    /**
     * @return the distance to the nearest region in this Regions object, or -1 
     *         if there is no region in the same contig, if a region onverlaps, 
     *         returns 0
     */
    @CompileStatic
    int distanceTo(final Region r) { 
        return this.index.containsKey(r.chr) ? this.index[r.chr].distanceTo(r.range) : (int)-1
    }
     
    List<Object> getExtrasAtPosition(String chr, int position) {
        List result = []
        
        RangeIndex chrIndex = this.index[chr]
        if(chrIndex == null)
            return result
        for(groovy.lang.Range range in chrIndex.getOverlaps(position)) {
            if(range instanceof GRange)
                result << range.extra
        }
        return result
    }
    

    /**
     * Add the specified range to this Regions file.
     * <p>
     * NOTE: the 'end' is treated as exclusive of the range covered. This
     * is consistent with BED file notation, but different to Groovy ranges.
     * Since most methods in this class use Groovy ranges, you will generally
     * get a Range object out that has end one less than the end you put in
     * with add. Thus if you are iterating through one BED file and adding the
     * ranges to another with this method, you <i>must</i> add one to the
     * end position that you pass to this method!
     * <p>
     * Will trigger the flag to indicate this BED is loaded in memory.
     */
    @CompileStatic
    Regions addRegion(String chr, int start, int end, Object extra = null) {
        
        // Add to full range index
        RangeIndex chrIndex = index[chr]
        if(!chrIndex) {
            chrIndex = new RangeIndex()
            index[chr] = chrIndex
        }
        
        groovy.lang.IntRange newRange = null
        newRange = new GRange(start,end-1,extra)
        chrIndex.add(newRange)
        
        if(!allRanges.containsKey(chr))
            allRanges[chr] = new ArrayList(1000)
            
        allRanges[chr] << newRange
        return this
    }
    
    /**
     * Adds the region to this Regions object in such a way that the same Region
     * object is returned in iteration, preserving any expando properties set on
     * the region.
     * 
     * @param r
     * @return
     */
    @CompileStatic
    Regions addRegion(Region r) {
        // Add to full range index
        String chr = r.chr
        RangeIndex chrIndex = index[chr]
        if(!chrIndex) {
            chrIndex = new RangeIndex()
            index[chr] = chrIndex
        }
        
        IntRange range = null
        if(r.range instanceof GRange && ((GRange)r.range).extra == null) {
            ((GRange)r.range).extra = r
            range = r.range
        }
        else {
            range = new GRange((int)r.range.from, (int)r.range.to, (Object)r)
        }
        
        chrIndex.add(range)
        
        if(!allRanges.containsKey(chr))
            allRanges[chr] = new ArrayList(1000)
        allRanges[chr] << range
        return this
    }
        
    /**
     * Simplify all overlapping regions down to a single region, with an optional
     * closure as a callback to combine the attributes of combined regions
     * <p>
     * The closure is passed two IntRange / GRange objects to combine. The 
     * callback shoud return the new {@link GRange#extra} object to 
     * assign to the result region.
     *
     * @return a new Regions that consists of all the overlapping regions of
     *          this Regions merged together.
     */
    @CompileStatic
    Regions reduce(Closure reducer=null) {
        Regions result = new Regions()
        this.index.each { String chr, RangeIndex ranges ->
            result.allRanges[chr] = []
            result.index[chr] = ranges.reduce(reducer)
            result.index[chr].each { IntRange r ->
                result.allRanges[chr].add(r)
            }
        }
        return result
    }

    /**
     * Return true iff the given chromosome and position fall into
     * a range covered by this BED file
     */
    @CompileStatic
    boolean contains(String chr, int position) {
        RangeIndex chrIndex = index[chr]
        return chrIndex != null && (position in chrIndex)
    }
    
    /**
     * Return a new Regions that contains the same regions as this one,
     * but ensuring each identical region is only present once.
     *
     * @return
     */
    Regions unique() {
        Regions result = new Regions()
        this.eachRange(unique:true) { chr, from, to, extra ->
            result.addRegion(chr, from, to+1, extra)
        }
        return result
    }
    
    boolean isEmpty() {
        return allRanges.isEmpty()
    }

    @Override
    @CompileStatic
    public Iterator<Region> iterator() {
        return new RegionIterator(this)
    }
    
    @CompileStatic
    Regions subtract(Regions other) {
        final Regions result = new Regions()
        allRanges.each { final String chr, final List<IntRange> ranges ->
            for(final IntRange range : ranges) {
                final List<IntRange> subtracted = other.subtractFrom(chr, range.from, range.to)
                for(IntRange r : subtracted) {
                    result.addRegion(chr, r.from, r.to+1)
                }
            }
        }
        return result
    }
    
    @CompileStatic
    Object getAtIndex(Object obj) {
        if(obj instanceof Integer) {
            int index = (int)obj
            int total = 0
            Map.Entry<String,List<IntRange>> entry = allRanges.find {
                total += it.value.size()
                total >= index
            }
            
            Range result = entry.value[index-total]
            if(result instanceof GRange && ((GRange)result).extra instanceof Region) {
                return ((GRange)result).extra
            }
            else {
                return new Region(entry.key, result)
            }
        }
        else
        if(obj instanceof List) {
            List<Region> result = []
            List<Integer> indices = (List<Integer>)obj
            int total = 0
            List offsets = [0] + allRanges.collect { rangeEntry -> total += rangeEntry.value.size();  }
            List keys = allRanges.collect { it.key }
            for(int index in indices) {
                def offsetIndex = offsets.findLastIndexOf { Integer offset -> offset <= index }
                def key = keys[offsetIndex]
                Range r = allRanges[key][index-offsets[offsetIndex]]
                if(r instanceof GRange && ((GRange)r).extra instanceof Region)
                    result.add((Region)((GRange)r).extra)
                else
                    result.add(new Region(key, r))
            }
            return result
        } 
    }
    
    @CompileStatic
    Object getAt(Object obj) {
        this.getAtIndex(obj)
    }
    
    int getNumberOfRanges() {
       def result = allRanges.collect { it }.sum { it.value.size() }
       if(result == null)
           return 0
       else
           return result
    }
    
    @CompileStatic
    Regions plus(Regions other) {
        Regions result = new Regions((Iterable<IRegion>)this)
        for(Region r in other) {
            result.addRegion(r)
        }
        return result
    }
    
    String toString() {
        int numRanges = getNumberOfRanges()
        if(numRanges>0) {
            "${numRanges} regions starting at ${this.iterator().next()}"
        }
        else {
            "Empty Region Set (no regions)"
        }
    }
    
    void save(String fileName) {
        save(null, fileName)
    }
    
    /**
     * Save the regions in BED format. If an 'extra' option is provided, this is called
     * as a closure to return an id field to use for each region.
     * 
     * @param options
     * @param fileName
     */
    @CompileStatic
    void save(Map options, String fileName) {
        Utils.withWriters([fileName]) { Writer w ->
            save(options,w)
        }
    }
    
    /**
     * Save the regions in BED format. If an 'extra' option is provided, this is called
     * as a closure to return an id field to use for each region.
     * 
     * @param options   extra: closure returning data to write into 4th column, 
     *                  sorted: true to sort lexically, Comparator to sort custom
     *                          (see {@link NumericRegionComparator)
     * @param fileName
     */ 
    void save(Map options, Writer w) {
        
        def regionsToSave = this
        if(options?.sorted) {
            
            if(options.sorted instanceof Comparator)
                regionsToSave = this.toSorted(options.sorted)    
            else
                regionsToSave = this.toSorted(new RegionComparator())    
        }
        
        Closure c = options?.extra
        
        if(c != null) {
            regionsToSave.each { w.println([it.chr, it.from, it.to+1, c(it)].join('\t')) }
        }
        else {
            regionsToSave.each { w.println([it.chr, it.from, it.to+1].join('\t')) }
        }
    }
    
    /**
     * Return a new regions object that has each distinct region of this object
     * with a "coverage" value assigned.
     * 
     * @return
     */
    Regions coverage() {
        Regions result = new Regions()
        this.index.each { String chr, RangeIndex ranges ->
            for(GRange r in ranges.coverage()) {
                result.addRegion(new Region(chr, r))
            }
        }
        return result        
    }
    
    List<Map> bkr() {
        this.collect { [ chr: it.chr, from: it.from, to: it.to] }
    }
    
    @CompileStatic
    Regions enhance() {
        Regions result = new Regions()
        for(Region r in this) {
            if(r instanceof GRange) {
                GRange gr = ((GRange)r.range)
                gr.extra = (Object)r;
                result.addRegion(r)
            }
            else {
                GRange gr = new GRange(r.from, r.to, null)
                Region newRegion = new Region(r.chr, gr)
                gr.extra = newRegion
                result.addRegion(newRegion)
            }
        }
        return result
    }
    
    /**
     * A convenience method to return the contained regions as a list of Map objects.
     * <p>
     * This is primarily aimed at interactive use in shell or notebook environments,
     * it is not performant or memory efficient.
     */
    @CompileStatic
    List<Map<String,Object>> toListMap() {
        (List<Map>)this.collect { Region r ->
            [
                chr: r.chr,
                start: r.from,
                end: r.to,
                span: r.size()
            ] + r.properties.grep { Map.Entry e -> (e.value instanceof String) || (e.value instanceof Number) }.collectEntries()
        }
    }
    
    /**
     * Select the given number of ranges from these, approximately evenly spaced
     * 
     * @param desiredRanges target number of ranges to preserve
     * @param minRanges     the minimum number of ranges to preserve on each chromosome
     * @return
     */
    @CompileStatic
    Regions thin(int desiredRanges, int minRangesPerChromosme) {
        int currentTotalRanges = this.numberOfRanges
        
        double proportionRetained = ((double)desiredRanges) / currentTotalRanges
        if(proportionRetained == 0) 
            throw new IllegalArgumentException("Proportion of regions retained is too small: $proportionRetained ($desiredRanges / $currentTotalRanges")
            
        int keepOneEvery = (int)(1 / proportionRetained)
        int i=0
        
        return this.grep { Region r -> 
            
            if(this.index[r.chr].numRanges*proportionRetained < minRangesPerChromosme) {
                return true
            }
            
            return (++i) % keepOneEvery == 0 
        } as Regions
    }
    
    /**
     * Returns the total span from the beginning of the first region to the end of 
     * the last region on the given contig (chromosome).
     * 
     * @return
     */
    @CompileStatic
    Region getSpan(String contig) {
        RangeIndex index = this.index.get(contig)
        if(index.is(null)) 
            return null;
        
        return new Region(contig, index.first().from..index.last().to)
    }

    @CompileStatic
    Regions getContigRegions(final String chr) {
        Regions result = new Regions()
        if(!this.index.containsKey(chr))
            return result

        for(IntRange r in this.index[chr]) {
            result.addRegion((Region)((GRange)r).extra)
        }
        return result
    }
    
    /**
     * Divide these regions into two regions objects with approximately the same
     * total bp in each half
     *
     * @return
     */
    @CompileStatic
    List<Regions> balancedSplit() {
        
        final long totalBp = this.size()
        final int halfBp = (int)(totalBp >> 1)
        final Regions firstHalf = new Regions()
        Regions secondHalf = new Regions()

        long accumulatedBp = 0
        Iterator<Region> iter = this.iterator()
        while(accumulatedBp<halfBp) {
            Region r = iter.next()
            accumulatedBp += r.size()
            
            if(accumulatedBp<halfBp)
                firstHalf.addRegion(r)
            else
                secondHalf.addRegion(r)
        }
        
        while(iter.hasNext()) {
            secondHalf.addRegion(iter.next())
        }
        
        return [firstHalf, secondHalf]
    }
}
