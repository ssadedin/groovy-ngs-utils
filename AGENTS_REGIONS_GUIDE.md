# Genomic Regions API Reference (LLM Guide)

Package: `gngs`. Groovy (2.x syntax). Prefer `@CompileStatic` everywhere possible.

## Core Types

### IRegion (interface)
- `String getChr()` — contig/chromosome name
- `IntRange getRange()` — positional range (inclusive both ends)

### Region extends Expando implements IRegion, Serializable
Heavyweight region with arbitrary expando properties. Inclusive of both endpoints.

**Construction:**
- `new Region("chr1:100-200")` — parse string
- `new Region("chr1", 100, 200)` — direct (both inclusive)
- `new Region(Map props, String chr, int from, int to)` — with metadata
- `new Region(IRegion other)` — copy
- `new Region(SAMRecord read)` — from alignment

**Properties:**
- `String chr`, `IntRange range`, `Integer from`, `Integer to`
- `Object extra` — attached data (from GRange.extra)
- Expando: set arbitrary props via `region.myProp = value`

**Methods:**
- `long size()` — number of positions spanned
- `boolean overlaps(IRegion other)` — any overlap?
- `boolean overlaps(String chr, int from, int to)` — any overlap?
- `boolean overlaps(Regions regions)` — overlaps anything in set?
- `boolean spans(IRegion r)` — fully contains other?
- `Region intersect(IRegion other)` — shared portion (returns EMPTY_REGION if none)
- `Region union(IRegion other)` — combined span (must overlap)
- `double mutualOverlap(IRegion other)` — min fractional overlap of either region
- `Region widen(int bp)` / `widen(int leftBp, int rightBp)` — expand boundaries (floor 0)
- `Region copy()` — shallow copy
- `boolean isCase(IRegion other)` — supports `r1 in r2` syntax (tests overlap, not containment)
- `int getMidpoint()` — midpoint position
- `boolean isEmpty()` — is EMPTY_REGION?
- `boolean isMinorContig()` — heuristic: alt haplotype / decoy / mito?
- `static boolean isMinorContig(String chr)` — static version
- `Region stripContigPrefix()` — remove "chr" prefix
- `String toString()` — `"chr1:100-200"`
- `String igv()` — IGV localhost goto URL

**Static:**
- `Region.EMPTY_REGION` — sentinel for empty results
- `static boolean overlaps(IRegion r1, IRegion r2)`
- `static Region union(IRegion r1, IRegion r2)`

### GRange extends IntRange implements Serializable
Lightweight range (no chr). Carries `Object extra` for metadata.

- `new GRange(int from, int to, Object extra)`
- `boolean spans(IntRange r)` — fully contains?
- `boolean overlaps(IntRange other)`
- `GRange intersect(IntRange other)` / `intersectRange(IntRange other)` — null if no overlap
- `GRange widen(int bp)` / `widen(int leftBp, int rightBp)`
- `static boolean overlaps(IntRange a, IntRange b)`
- `Object extra` — arbitrary attached data

### Regions implements Iterable\<Region\>
Collection of regions indexed by chr for fast overlap queries. BED-like semantics on add (end exclusive), but stored ranges are inclusive.

**Construction:**
- `new Regions()` — empty
- `new Regions(Iterable<IRegion> regions)` — from list
- `new Regions(String chr, Iterable<Range> ranges)`
- `[new Region("chr1:1-10"), new Region("chr1:5-15")] as Regions` — Groovy coercion

**Adding (end is EXCLUSIVE, BED convention):**
- `addRegion(String chr, int start, int end)` — end exclusive
- `addRegion(String chr, int start, int end, Object extra)`
- `addRegion(Region r)` — preserves expando properties

**Size/Count:**
- `long size()` — total bp across all regions (NOT count of regions)
- `int getNumberOfRanges()` — count of regions
- `boolean isEmpty()`

**Overlap queries:**
- `boolean overlaps(IRegion r)` — any overlap?
- `boolean overlaps(String chr, int from, int to)`
- `boolean overlaps(Regions other)` — any mutual overlap?
- `List<IntRange> getOverlaps(IRegion r)` — raw ranges overlapping
- `List<IntRange> getOverlaps(String chr, int start, int end)` — both ends inclusive
- `List<Region> getOverlapRegions(IRegion r)` — as Region objects
- `boolean isCase(IRegion r)` — supports `region in regions` syntax

**Set operations:**
- `Regions intersect(Regions other)` — pairwise intersections of all ranges
- `Regions intersectRegion(Region other)` — intersect all ranges with one region
- `List<IntRange> intersect(String chr, int start, int end)` — intersected ranges
- `Regions subtract(Regions other)` — remove other's ranges from this
- `List<IntRange> subtractFrom(String chr, int start, int end)` — remove this from interval (end exclusive)
- `List<Region> subtractFrom(Region region)` — remove this from region
- `Regions reduce(Closure reducer=null)` — flatten/merge overlapping ranges; optional closure combines extras
- `Regions plus(Regions other)` — combine two Regions

**Navigation:**
- `IntRange nextRange(String chr, int pos)` — next range after pos
- `IntRange previousRange(String chr, int pos)` — prior range before pos
- `IntRange nearest(String chr, int pos)` — closest range to pos
- `int distanceTo(String chr, int pos)` — bp to nearest range (-1 if no ranges on chr, 0 if overlapping)
- `int distanceTo(Region r)` — bp to nearest range
- `IntRange forward(String chr, int pos, int count=1)` — nth range forward
- `IntRange backward(String chr, int pos, int count=1)` — nth range backward
- `List<Region> window(Region r, int n)` — n regions upstream + downstream

**Positional queries:**
- `boolean contains(String chr, int position)` — position inside any range?
- `List<IntRange> startingAt(String chr, int pos)` — ranges starting exactly at pos
- `List<Region> regionsStartingAt(String chr, int pos)` — as Regions
- `List<Range> endingAt(String chr, int pos)` — ranges ending exactly at pos
- `List<Region> regionsEndingAt(String chr, int pos)` — as Regions

**Transformation:**
- `Regions widen(int bp)` — widen all regions
- `Regions uniquify()` — deduplicate identical ranges
- `Regions enhance()` — ensure all ranges have Region as extra
- `Regions thin(int desiredRanges, int minPerChr)` — evenly spaced subset
- `Regions coverage()` — breakpoint regions with overlap count as extra
- `Regions getContigRegions(String chr)` — subset for one chromosome
- `Region getSpan(String contig)` — min-to-max region on contig
- `List<Regions> balancedSplit()` — split into two ~equal bp halves

**Iteration:**
- Implements `Iterable<Region>` — use `.each{}`, `.collect{}`, `.grep{}`, `.find{}`, etc.
- `eachRange(Closure c)` — iterate as `(chr, from, to)` or `(chr, from, to, extra)` or `(chr, range)` or `(region)`
- `eachRange(Map options, Closure c)` — options: `unique:true`

**Indexed access:**
- `regions[0]` — nth region across all chromosomes
- `regions[1..3]` — sublist by index

**I/O:**
- `save(String fileName)` — write BED
- `save(Map options, String fileName)` — options: `extra:` closure for 4th col, `sorted: true` or Comparator

**Utility:**
- `double mutualOverlap(IRegion r)` — max mutual overlap of any contained region with r
- `List<Map> toListMap()` — convert to list of maps (interactive use)
- `List<Map> bkr()` — minimal list of maps `[chr, from, to]`
- `String toHTML(Map options=[:])` — HTML table with IGV links

### RangeIndex implements Iterable\<IntRange\>
Per-chromosome index. Usually accessed via `Regions.index[chr]`.

- `add(IntRange r)` / `add(int start, int end, Object extra=null)` — end exclusive on int version
- `remove(IntRange r)` — remove exact range
- `List<IntRange> getOverlaps(int pos)` — ranges at position
- `List<IntRange> getOverlaps(int start, int end)` — overlapping ranges (inclusive)
- `boolean overlaps(int start, int end)` — any overlap?
- `List<IntRange> intersect(int start, int end)` — clipped to query
- `List<IntRange> subtractFrom(int start, int end)` — gaps after removing index ranges
- `List<IntRange> startingAt(int pos)` / `endingAt(int pos)`
- `IntRange nearest(int pos)` / `int distanceTo(int pos)`
- `Range nextRange(int pos)` / `previousRange(int pos)`
- `IRangeIndex reduce(Closure reducer=null)` — merge overlapping
- `Iterator<IntRange> iterator()` / `iteratorAt(int pos)` / `reverseIterator()` / `reverseIteratorAt(int pos)`
- `int getNumRanges()`

## Comparators

- `RegionComparator` — lexical chr, then positional
- `NumericRegionComparator` — numeric chr (handles chrX/chrY), then positional

## Key Patterns

```groovy
// Create regions
regions = [new Region("chr1:100-200"), new Region("chr1:150-300")] as Regions

// Programmatic add (end EXCLUSIVE)
regions.addRegion("chr1", 100, 201)

// Overlap check
if(region in regions) { ... }
if(regions.overlaps(region)) { ... }

// Get overlapping regions
List<Region> hits = regions.getOverlapRegions(myRegion)

// Flatten overlaps
Regions flat = regions.reduce()

// Subtract
Regions diff = regionsA.subtract(regionsB)

// Intersect
Regions shared = regionsA.intersect(regionsB)

// Iterate with Groovy power
regions.grep { it.chr == "chr1" && it.size() > 1000 }
regions.collect { [it.chr, it.from, it.to] }
regions.sum { it.size() }
```

## Gotchas

1. `Regions.addRegion(chr, start, end)` treats end as **exclusive** (BED convention), but stored IntRange is **inclusive** (end-1).
2. `Region("chr1:100-200")` — both 100 and 200 are **inclusive**.
3. `regions.size()` returns **total bp**, not count. Use `regions.numberOfRanges` for count.
4. `Region.EMPTY_REGION` is a sentinel — test with `.isEmpty()` or `.is(EMPTY_REGION)`.
5. GRange.extra often holds the parent Region — access via `((GRange)range).extra`.
6. Use `@CompileStatic` and avoid Groovy 3.x Parrot syntax (must compile under 2.x).
