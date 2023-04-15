# Introduction

Sequence alignments are central to working with genomic data. Groovy NGS supplies a wrapper class
that adds idiomatic Groovy constructs to the underlying HTSJDK and Picard classes that provide
Java support for working with genomic data. Unlike the support for genomic regions, this layer
is a relatively thin wrapper since HTSJDK already provides a very effective API.

# Opening BAM and CRAM Files

The core class for interacting with alignment files is the [SAM](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/SAM.html)
class. This class *only* works with indexed BAM, CRAM or SAM files. It provides a wealth of functions that allow
highly efficient, streamlined access to reads within these formats.

To creater a `SAM` object, just use the constructor:

```groovy
bam = new SAM("test.bam")
```

To see which contigs are in the file, you can access the `contigs` property:

```````columns
left:
**Input**
```groovy
bam.contigs
```

right:
**Output**
`````text
chrM=16571
chr1=249250621
chr2=243199373
chr3=198022430
...
`````
````````

To see which samples are in an alignment file, use the `samples` property:

```````columns
left:
**Input**
```groovy
bam.samples
```

right:
**Output**
`````text
[NA12878, NA12878]
`````
````````

Note that each sample may be returned multiple times if there is more than one read group for the sample
in the alignment file.


# Accessing Reads

Reads within the `SAM` file are accessed using several different methods. The right method depends on

- whether you want to access reads from a specific region
- if you need access to both pairs of paired end reads at the same time

For basic access, use the [eachRecord](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/SAM.html#eachRecord(java.util.Map,%20groovy.lang.Closure))
method to iterate over every read in the SAM file, one at a time:

```````columns
left:
**Input**
```groovy
count = 0
bam.eachRecord { ++count }
println "There are $count reads in the file"
```

right:
**Output**
`````text
There are 1899 reads in the file
`````
````````

Within the closure, a [SAMRecord](https://www.javadoc.io/doc/com.github.samtools/htsjdk/1.133/htsjdk/samtools/SAMRecord.html) object is passed
as an argument and all of its properties are accessible.

For example, the maximum insert size could be calculated using the `inferredInsertSize` property of each `SAMRecord` object:

```````columns
left:
**Input**
```groovy
maxInsertSize = 0
bam.eachRecord { r ->
    maxInsertSize = Math.max(maxInsertSize, Math.abs(r.inferredInsertSize))
}
println "The maximum insert size is $maxInsertSize"
```

right:
**Output**
`````text
The maximum insert size is 114074
`````
````````

While `eachRecord` is a simple method, often a more idiomatic and functional style is achievable
by using the [withIterator](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/SAM.html#withIterator(groovy.lang.Closure)) 
which makes available a Groovy [Iterator](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Iterator.html), affording
all of the usual Groovy enhancements that are available on iterator objects.

For example, the same can be achieved with an iterator like so:

```````columns
left:
**Input**
```groovy
bam.withIterator { it*.inferredInsertSize.max() }
```
right:
**Output**
`````text
114074
`````
````````

## Accessing reads from a specific locus

To access reads from a given locus, you can supply a `Region` object as a parameter to the `withIterator` or `eachRecord` methods to 
specify the region.


# Accessing Paired Reads

It is often useful to access reads in their pairs. However, within coordinate sorted alignment files, paired reads
may be spaced far apart which makes retrieval of the mate for each read potentially very expensive. Groovy NGS offers
an optimised method, [eachPair](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/SAM.html#eachPair(groovy.lang.Closure)).

```````columns
left:
**Input**
```groovy
count = 0
bam.eachPair { r1, r2 -> ++count}
println "There are $count read pairs"
```
right:
**Output**
`````text
There are 942 read pairs
`````
````````

Note: accessing read pairs can be memory intensive, as each read is buffered until its mate is found. Unpaired reads
will not be output.

# Generating Pileups

Genomic data often contains many reads covering any one position and it is of interest to 
access the full set of reads covering the position (known as a "pileup"). Groovy NGS offers several
methods on the `SAM` class as well as a dedicated [Pileup](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/Pileup.html)
class to make generating pileups easier.

## Pileup at a Specific Position

To generate a pileup at a specific position, use the [pileup](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/SAM.html#pileup(java.lang.String,%20int))
method. This returns a `Pileup` object that can then be queried for summary information as well as the full list of reads covering the position:

```````columns
left:
**Input**
```groovy
p = bam.pileup('chr1',16890560)
println "There are ${p.alignments*.read.size()} reads at the position"
println "There are ${p.summaryAsMap.C} C bases at the position"
```
right:
**Output**
`````text
There are 3 reads at the position
There are 3 C bases at the position
`````
````````

From this we can see there are 3 reads which all have the base `C` at that position.

## Pileups at Multiple Positions

Generating a pileup at a single position is very compute intensive as it requires fully scanning
the context around the position to find all the reads. If you are analysing a whole region, you 
should instead create a `PileupIterator` over the whole region using the [pileup(chr,start,end)](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/SAM.html#pileup(java.lang.String,%20int,%20int))
method. This will return an iterator that you can apply standard Groovy methods to enumerate the
pileup state efficiently at each position within the range.

```````columns
left:
**Input**
```groovy
bam.pileup('chr1',16890550, 16890560 ).collect { p ->  p.alignments.size() }
```
right:
**Output**
`````text
[2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3]
`````
````````

# Calculating Coverage Depth

As seen above, coverage can be calculated using the `pileup` method on the `SAM` class. There are several
direct methods for calculating coverage depth which are more efficient, however, including some that
directly compute statistics for you:

```````columns
left:
**Input**
```groovy
bam.coverageStatistics('chr1',16890550, 16890560 ).mean
```
right:
**Output**
`````text
2.090909090909091
`````
````````

Note that the object returned is an Apache commons-math [SummaryStatistics](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/stat/descriptive/SummaryStatistics.html)
object, allowing you to access any of its properties.
