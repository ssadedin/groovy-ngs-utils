# Regions and Ranges

DNA is most commonly conceptualised as a set chromosomes or contigs represented
by linear coordinate systems. Defining and manipulating spans within these
spaces is therefore been critical to working with genomic data. Entities such as
genes, exons, transcripts, variants, and many similar constructs are usually
represented as coordinates within the genome, typically referred to as spans,
ranges or regions. This section explains how to work with genomic regions in
Groovy-NGS.

Groovy itself provides interfaces and a small amount of
native support for working with ranges of values, such as the [Range](https://docs.groovy-lang.org/latest/html/api/groovy/lang/Range.htm) class. Unfortunately this support is very minimal, so Groovy-NGS builds on these to provide full featured range operation support. In general, the terminology is used as follows:

## Terminology

- **Range** : A contiguous sequence of genomic positions not qualfied by a contig / chromosome (defined by start and end position)
- **Region** : A range that is further associated with a specific sequence, ie: a chromosome, or contig. Operations treat Regions associated with different sequences as independent entities that do not interact.
               
Note: for historical reasons, Groovy-NGS uses `chr` to refer to the contig or chromosome with which a 
`Region` is associated. For all intents and purposes, you may think of the `chr` as a contig, or even
an abstract sequence name.

## Key Supporting Classes 

Formally, the key classes supporting range operations in Groovy-NGS are:

```table {title: "Range Data Types"}
Type, Description
[IRegion](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/IRegion.html), Core interface implemented by all data types that can interface with range operations - supports only contig and interval. Data types representing entities that are localisable to a specific position in the genome usually implement `IRegion` to specify their location.
[Region](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/Region.html), "A region implementation that is also a groovy [expando](https://docs.groovy-lang.org/latest/html/api/groovy/util/Expando.html). This means that you can set your own properties directly on this object, allowing convenient association of meta data or ancilliary information to the region in data analysis."
[Regions](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/Regions.html), "A collection of potentially overlapping `Region` objects, stored in indexed fashion for rapid lookup."
```

## Implementation

Under the hood, Groovy-NGS stores each range in a `Regions` object using an
index that enables efficient lookup. Unlike the commonly used approach of applying
an [Interval Tree](https://en.wikipedia.org/wiki/Interval_tree), Groovy-NGS uses
a boundary index. That is, when a new range is inserted, the two end points are inserted
into a tree map ordered on position, specific to the chromosome / contig. At places where
the range overlaps existing ranges already in the index, breakpoint entries are already
inserted. This functionality is supported by the [RangeIndex](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/RangeIndex.html) class.

Note: the breakpoint index design allows for extremely rapid local lookup, but
creates
some worst case performance outcomes you may wish to look out for. In particular,
the scenario where many small regions are inserted first and then a set of large regions
are inserted after causes the insertion operation to become O(N) on the number
of small regions. Now that HTSJDK offers a native
[IntervalTree](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/util/IntervalTree.html)
an implementation backed by this class may be considered to allow this worst
case to be avoided. If you encounter the worst case described here, using the
HTSJDK `IntervalTree` class instead may be an option.


