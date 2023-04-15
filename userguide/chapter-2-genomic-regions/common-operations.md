
This section describes how to do a variety of common operations that are useful in
analysing genomic data.

# Creating Regions Objects

The first thing you need to do to work with regions is to create a [Regions](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/Regions.html) object.
One way to do this is "by hand" - ie: by manually constructing [Region](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/Region.html) objects.
You can do this by passing a region defined as a string using the standard format to
the `Region` constructor, for example:

```groovy
r = new Region('chr1:100-200')
```

Then you can convert a list of these to a `Regions` object using the Groovy `as` 
operator. For example, consider the following set of regions:

```
          A
|--------------------|
      B       C            D
    |---|   |---|        |---|
0   5   10  15  20   25  30  35
```

To create a `Regions` object representing these directly, we can use the following code:

```groovy
regions = [
              'chr1:5-10', 'chr1:15-20', 'chr1:30-35', 'chr1:0-25'
          ].collect { 
              new Region(it) 
          } as Regions
```

Observe that in this set of example regions we have overlapping ranges. In
general, Groovy-NGS does not automatically combine or merge regions, and the resulting
`Regions` object behaves much like a simple collection of the regions you put in, only
with some operations accelerated. If you desire that the regions be combined, you 
need to use one of the specific operations to do that after creating your `Regions` object.

# Basic Metrics

Some simple metrics you might like to know about your set of include:

| Description | Code   | Output |
|-|-|-|
| the total span of all the ranges included | `regions.size()` |  `44` |
| the count of the number of ranges included | `regions.numberOfRanges` |  `4` |

Note: the use of `size()` slightly conflicts with the standard Groovy semantics 
for the `size()` method as this would normally return the number of elements in a collection. Instead, you can think of this as returning the total number of bases
in the collection.


# Treating as a Collection

The `Regions` object implements standard Java and Groovy interfaces and methods to enable
common collections operations. In particular, the [Iterable<Region>](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/lang/Iterable.html) interface is implemented so that the full range of Groovy 
enhanced operations for collections are available. These operations are extremely powerful and
often allow you to accomplish complex data analysis tasks with just a few lines of chained
declarative invocations. When used as an `Iterable`, the `Regions` class will
span across all chromosomes / contigs. If you want to operate at the
chromosome/contig level, you can consider using the
[index](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/Regions.html#index)
attribute to access the per-chromosome index which also supports the Iterable interface.

## Indexed Access

For accessing subsets of regions within a `Regions` object, specific support is provided to
allow lookup using "square bracket" notation efficient. Therefore, for example, to access
the 2nd and 3rd regions of our set, the following code will work:

**Input**
```groovy
regions[1..2]
```

**Output**
```
[chr1:15-20, chr1:30-35]
```


# Finding Overlaps

One of the most common needs is to find which regions overlap a given interval. Generally,
the most practical method to use for this is the [getOverlapRegions](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/Regions.html#getOverlapRegions(gngs.IRegion)) method. This method will return a list of the regions from the set that have at least
1 base of overlap with the region that is specified.

Example: to find which regions overlap the range 12 to 30:

**Input**:

```groovy
r = new Region('chr1:12-30')
overlaps = regions.getOverlapRegions(r)
```

**Output**:

```
[chr1:0-25, chr1:15-20, chr1:30-35]
```

Note: the `getOverlapRegions` method returns a `List` object rather than a `Regions` object.
Therefore in this case, the `.size()` method returns the number of elements, not the total
span of all the regions contains. In general, you will find Groovy NGS returns 
plain collections as a result of its methods, which you can then choose to convert into
a `Regions` object by adding `as Regions` if that is what you need.

Groovy NGS interprets region boundaries as *inclusive*. This means that the region
`chr1:0-10` includes the position 10 and does overlap the region `chr1:10-20`.

Note: There are several other `getOverlaps` methods. These can be slightly more efficient, in specific circumstances.

# Intersection

Intersection finds the shared parts of each region in your set with another region that 
you specify. To this, use the [intersectRegion](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/Regions.html#intersectRegion(gngs.Region)) method. 

**Input**
```groovy
ix = regions.intersectRegion(r)
ix*.toString()
```

**Output**

```groovy
[chr1:12-25, chr1:15-20, chr1:30-30]
```

# Flattening

In many types of analysis, you desire to know only the total set of regions encompassed
by an overlapping set. Groovy-NGS mirrors the R-lang terminology here, offering a `reduce`
function to flatten a set of regions to an equivalent set where all overlaps are compressed
to a single region.

**Input**

```groovy
regions.reduce()*.toString()
```

**Output**

```
[chr1:0-25, chr1:30-35]
```

This has changed the ranges to the following:

```
|--------------------|   |---|
0   5   10  15  20   25  30  35
```



# Assigning Properties

In analysing genomic data we usually are interested in specific attributes of the regions we are
working with. For example, the genes they represent, statistics about their attributes or any arbitrary
data that we associate with the region for the purpose of the analysis. It is very convenient therefore to be able to attach metadata to each region so that these are available at any point
they are needed in analysis.

The Groovy-NGS
[Region](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/Region.html) class
extends the Groovy
[Expando](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/Region.html)
class which allows you to set arbitrary properties on the object alongside the properties
that the class natively supports. For example:

**Input**

```groovy
regions[0].foo = 'bar'
regions[1].foo = 'baz'
regions[2].foo = 'boo'
regions[3].foo = 'boz'

regions.each { r ->
    println "The region $r has foo = $r.foo"
}
```

**Output**

```
The region chr1:5-10 has foo = bar
The region chr1:15-20 has foo = baz
The region chr1:30-35 has foo = boo
The region chr1:0-25 has foo = boz
```

There are some important caveats with using this feature:

- you should ensure the properties assigned don't clash with native property names. In such cases,
  which property is accessed may be ambiguous.
- some operations on `Regions` objects that modify the regions will lose the metadata
  because they return new `Regions` objects. If you want to re-associate the metadata
  you may have to do this manually (for example, by using `findOverlapRegions` or similar.
  

