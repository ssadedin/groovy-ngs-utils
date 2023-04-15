# Loading BED files

Regions are commonly loaded from and saved to BED files. Groovy NGS provides the 
[BED](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/BED.html) class to read data 
from BED formatted files.

The BED class is designed to work with both in-memory BED files and also provide a subset of 
operations that can work in streaming fashion with files on storage. Therefore its constructor
accepts a file path but will not load the contents automatically. To load a `Regions` object into 
memory from a BED file, use the 
[load](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/BED.html#load(java.util.Map))
method. Consider this bed file:

```text {title: "test.bed"}
chr1	100	150  r1
chr2	190	250  r2
chr3	300	350  r3
```

Then we can load it to a regions object with:

```groovy
bed = new BED("test.bed").load()
assert bed.numberOfRanges == 3
```

This will load a minimalist version of the BED that does not include the
fourth column (usually referred to as the `name` or `id` column in the BED file).
This  is designed to minimise memory use and the resultant `Regions` object 
does not support the assignment of additional properties. To load a 
fully functional `Regions` object, add the `full:true` optional parameter
in either the constructor or the `load` function:

```groovy
bed = new BED("test.bed", full:true).load()
assert bed.numberOfRanges == 3
assert bed[0].id == "r1"
```

# Saving BED files

To save data in BED format from a `Regions` object, use the [Regions.save](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/Regions.html#save(java.lang.String)) method on the `Regions` object.

By default, this will save the regions in the order they were loaded, without any additional
columns. To add the `id`/`name` column as the fourth column, specify a closure that returns the
value to store there using the `extra` attribute. For example, to save the above BED file
with the same contents to a new file name, you would use:

```groovy
bed.save('test2.bed', extra: { it.id })
```

It's common to want to save regions in sorted order (eg: so that the file can be indexed). To do this,
supply the `sorted` attribute. If you supply the value `true` then the file will be written
in sorted order using lexical sorting on the contig/chromosome name. This means that `chr10` will follow
directly after `chr1`.  You can sort in
a custom order by setting the `sorted` attribute to your own implementation of a `Comparator<Region>`.

If you would like to store chromosomes in genomic order (ie: `chr2` follows `chr1`), Groovy NGS provides
the [NumericRegionComparator](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/NumericRegionComparator.html) class
that you can use. For example:

```groovy
bed.save('test2.bed', extra: { it.id }, sorted: new NumericRegionComparator())
```

# Region Based Data Tables

It is common to encounter files stored in column based format (for example,
tab separated) with genomic coordinates in some of the columns and other data
in the remaining columns.

To load this kind of data, use the [RangedData](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/RangedData.html)
class. This class expects a tab or comma separated file with headers defining the column names. The class
will auto-detect data types for columns based on the contents of the first rows (numeric columns being 
parsed as numbers). By default the first 3 columns are interpreted as the contig / chromosome, start 
coordinate and end coordinate repsectively. However you can change this by providing the column indices 
for these in the constructor.

Consider the following file:

```text {title: "test.tsv"}
chr	start	end	name	age
chr1	100	120	simon	10
chr1	140	210	fred	15
chr1	190	250	jane	20
chr1	300	350	tom	13
```

We can load a `Regions` object with `name` and `age` populated as follows:

```
regions = new RangedData("test.tsv").load()

assert r[0].name == "simon"
assert r[0].age == 10
```

For more details, see the [RangedData](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/RangedData.html) class
API documentation.

