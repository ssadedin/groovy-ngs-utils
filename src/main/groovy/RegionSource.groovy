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

import groovy.transform.CompileStatic;

class RegionSource implements Iterable<Region> {
		
	/**
	 * Index for looking up overlaps
	 */
	Map<String, RangeIndex> index = [:]
	
	/**
	 * A list of ranges in the order they were loaded
	 */
	Map<String,List> allRanges = [:]
	
	/**
	 * Map containing all the ranges loaded from the bed file,
	 * index by start position, and partially reduced to make
	 * lookups efficient.
	 * Note: only loaded if load() is called!
	 */
	TreeMap<Long,Range> startRanges = new TreeMap()
	
	/**
	 * Map containing all the ranges loaded from the bed file,
	 * indexed by end position, and partially reduced to make
	 * lookups efficient.
	 */
	TreeMap<Long,Range> endRanges = new TreeMap()
	
	/**
	 * Empty bed file
	 */
	RegionSource(Map attributes=[:]) {
	}
	
	/**
	 * Create a BED from a list of regions
	 */
	RegionSource(Map attributes=[:], Iterable<Region> regions) {
		for(Region r in regions) {
			addRegion(r.chr, r.from, r.to)
		}
	}
	
	@CompileStatic
	int size() {
		int sizeCount = 0
		this.eachRange { String chr, int start, int end ->
			sizeCount += (end - start + 1)
		}
		return sizeCount
	}
	
	void eachRange(Closure c) {
	   eachRange(unique:false, c)
	}
	
	/**
	 * Iterate over each line in the BED file and invoke the specified
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
			  else {
				  c.call(chr,r.from,r.to)
			  }
				  
			  progress.count()
			}
		}
	}
	
	List<Range> intersect(String chr, int start, int end) {
		RangeIndex chrIndex = this.index[chr]
		if(chrIndex == null)
			return []
		
		return chrIndex.intersect(start,end)
	}
	
	List<Range> getOverlaps(String chr, int start, int end) {
		RangeIndex chrIndex = this.index[chr]
		if(chrIndex == null)
			return []
		return chrIndex.getOverlaps(start,end)
	}
	
	List<Range> subtractFrom(String chr, int start, int end) {
		RangeIndex chrIndex = this.index[chr]
		if(chrIndex == null)
			return [start..end-1]
		
		return chrIndex.subtractFrom(start,end)
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
	
	List<Range> startingAt(String chr, int pos) {
		RangeIndex chrIndex = this.index[chr]
		if(chrIndex == null)
			return []
		
		return chrIndex.startingAt(pos)
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
	groovy.lang.IntRange nextRange(String chr, int pos) {
		RangeIndex chrIndex = this.index[chr]
		if(chrIndex == null)
			return null
		return chrIndex.nextRange(pos)
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
	 * Add the specified range to this RegionSource file.
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
	RegionSource addRegion(String chr, int start, int end, Object extra = null) {
		
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
	 * Simplify all overlapping regions down to a single region
	 *
	 * @return a new BED that consists of all the overlapping regions of
	 *          this BED merged together.
	 */
	RegionSource reduce() {
		RegionSource result = new RegionSource()
		eachRange { String chr, IntRange r ->
			
			// Already processed this region
			if(result.getOverlaps(chr, r.from, r.to))
				return
				
			// Find the overlaps
			List<Range> overlaps = getOverlaps(chr, r.from, r.to)
			IntRange mergedRange = new IntRange(overlaps.min { it.from }.from, overlaps.max { it.to }.to)
			result.addRegion(chr, mergedRange.from, mergedRange.to+1)
		}
		return result
	}

	/**
	 * Return true iff the given chromosome and position fall into
	 * a range covered by this BED file
	 */
	boolean contains(String chr, int position) {
		RangeIndex chrIndex = index[chr]
		return chrIndex != null && (position in chrIndex)
	}
	
	/**
	 * Return a new RegionSource that contains the same regions as this one,
	 * but ensuring each identical region is only present once.
	 *
	 * @return
	 */
	RegionSource unique() {
		RegionSource result = new RegionSource()
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
		return new Iterator<Region>() {
			
			Iterator allIterator = RegionSource.this.allRanges.iterator()
			
			// Iterator chrIterator = allIterator.hasNext() ? allIterator.next().value.iterator() : null
			Iterator chrIterator =  null
			
			
			Map.Entry<String, RangeIndex> currentChr
			
			String chr
			
			boolean hasNext() {
				
				
				if(!chr)
					return RegionSource.this.allRanges.any { it.value.size() > 0 }
				return chrIterator.hasNext()
			}
			
			Region next() {
				if(chr == null)
					nextChr()
			   
				Range range =  chrIterator.next()
				if(!chrIterator.hasNext())
					nextChr()
				 
				return new Region(chr,range)
			}
			
			void nextChr() {
			   while(allIterator.hasNext()) {
				   currentChr = allIterator.next()
				   if(currentChr.value) {
					   chrIterator = currentChr.value.iterator()
					   chr = currentChr.key
                       break
				   }
			   }
			}
			
			void remove() {
				throw new UnsupportedOperationException()
			}
		}
	}
	
	int getNumberOfRanges() {
	   def result = allRanges.collect { it }.sum { it.value.size() }
	   if(result == null)
		   return 0
	   else
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
}
