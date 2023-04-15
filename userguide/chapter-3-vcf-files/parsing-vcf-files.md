Groovy NGS provides a dedicated VCF parser designed to offer highly streamlined access to VCF 
information, via the [gngs.VCF](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/VCF.html) class.
This parser allows access to key header properties as well as full information about variants.

# Loading a VCF File

To load a VCF into memory, use the [VCF.parse](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/VCF.html#parse(Closure)) 
static method. This method accepts a closure which, if provided, is called for each variant to determine if the variant
should be retained. This allows loading a subset of variants from a VCF. For example, here we load 
only variants on chromosome 22.

```````columns
left:
**Input**
```groovy
vcf = VCF.parse("test.vcf") {
   it.chr == 'chr22'
}
println "The samples are: $vcf.samples"
println "There are ${vcf.size()} variants on chr22"
```
right:
**Output**
`````text
The samples are: [NA12877, NA12878, NA12879]
There are 24 variants on chr22
`````
````````

# Accessing Header Information Only

If access to the variant information is not required, you can avoid parsing the whole VCF by calling
the constructor directly:

```````columns
left:
**Input**
```groovy
vcf = new VCF("test.vcf")
println "The samples are: $vcf.samples"
```
right:
**Output**
`````text
The samples are: [NA12877, NA12878, NA12879]
`````
````````

# Streaming Processing

For large VCFs it can be beneficial to process a VCF file line by line. In this case, the 
[VCF.filter](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/VCF.html#filter(java.lang.String,%20Closure)) method
can be used to process variants in streaming fashion. By default, the `filter` method reads from standard 
input and writes to standard output, but a file names can be provided to read from a file.

```groovy
VCF.filter('test.vcf') {
   it.chr == 'chr22' // output only variants on chr22 to stdout
}
```

# Processing Variants in VCFs

Each variant is represented by an instance of the [gngs.Variant](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/Variant.html) 
class which makes available detailed information about the variant. The `VCF` class  implements the `Iterable<Variant>` interface,
making available standard operations that are available on all Groovy [Iterators](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Iterator.html).

```````columns
left:
**Input**
```groovy
vcf = VCF.parse("test.vcf").countBy { it.type }
```
right:
**Output**
`````text
{SNP=22, DEL=1, INS=1}
`````
````````

# Indexed VCFs

When accessing a small number of variants from a large VCF, it is more efficient to index the VCF and parse
only the relevant lines. The [gngs.VCFIndex](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/VCFIndex.html) class
can be used to query variants within a specified genomic region.


```````columns
left:
**Input**
```groovy
index = new VCFIndex("test.vcf.gz")
index.iterator("chr1", 13792200, 13792300).collect { it }
```
right:
**Output**
`````text
[chr1:13792284 G/T]
`````
````````
