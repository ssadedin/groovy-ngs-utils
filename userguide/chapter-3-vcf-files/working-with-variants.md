# General Notes

Variants are one of the core focuses of many genomic analysis, so Groovy NGS goes to significant effort to
streamline and optimise access to variant information. This is one of the exceptions, where 
even though HTSJDK does provide a [VariantContext](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html) class,
Groovy NGS provides a different, substitute class as a complementary option. This is because  HTSJDK focuses on completely 
representing all features from the VCF specification and also ensuring high levels of safety and robustness in how it is used.
Groovy NGS, therefore, provides an alternative that is highly ergonomic, but makes some important simplifiying assumptions and 
only presents a commonly used subset of the VCF specification.

The most important simplifying assumption is that many operations, there is a "simple" version that returns results for the first 
allele and / or the first sample in a VCF. This assumption means that many analyses where these assumptions hold true 
can be expressed with vastly simpler code and especially for interactive analyses (such as in Jupyter Notebooks), this 
adds a great level of convenience.

Warning: you may miss variants when using default methods with multi-allelic variants. In general, you should
apply Groovy NGS classes to VCFs that have been "normalised" and reduced to primitives, for example, by using 
[bcftools norm](https://samtools.github.io/bcftools/bcftools.html#norm).

# Variant Properties

The key class for working with variants in Groovy NGS is the [Variant](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/Variant.html) 
class. Conceptually, this class represents a single line in a VCF file, and hence it can include multiple alleles 
and multiple samples at a single site.

## Alleles

The base sequence of default alternate allele is accessible as the `alt` property. This value is exactly
as presented in the first element of the ALT column from the VCF file. Therefore, for indels it includes
the context base that is generally present for those variant types. The full list of alternate allele 
bases is accessible via the `alts` property.

These properties tell you only the base sequence of the allele. For more information, you can access
a full list of `Allele` objects as the `alleles` property. These allow you to directly query
the type of the variant, start and end position and other attributes.

## Dosage

GNGS refers to the number of copies of a variant that are present as the `dosage`. Therefore, in a 
diploid genome, having genotype `1/1` will yield a dosage of `2` for the alternate allele, while
having genotype `1/0` will give dosage of 1. This generally makes calculations more straightforward
than dealing with the exact genotype, especially when accessing the dosage for the first sample and 
allele (directly as the `dosage` property).

To see the dosage of the default allele for every sample, you can query the `dosages` property.
For example assuming the proband is the first sample in a trio VCF, we can identify de novo variants
as follows:

```groovy
vcf.grep {
    it.dosages[0] > 0 && it.dosages[1]==0 && it.dosages[2] == 0
}
```

## Depth

The total read depth at the site of a variant can be queried using the `totalDepth` property. To
get the depth of the first allele, you may query the `altDepth` property. To see the alternate 
allele depths for other alleles,  you can call `getAlleleDepths(alleleIndex)` where providing
0 specifies the reference allele.

Often we are interested in the proportion of reads supporting a variant (refererred to as 
the `variant allele frequency`. This is returned as the `vaf` property for the default 
allele and sample. To query it for arbitrary alelel and sample, call the
[getVaf(alleleIndex, sampleIndex)](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/Variant.html#getVaf(int,%20int))
method.

# Variants as Regions

Like other entities that belong to a specific location in the genome, Variants implement the
`IRegion` interface which allows variants to be passed to many functions for region based
operations as a parameter.

For example, to test if a particular variant lies within a region of interest:

```````columns
left:
**Input**
```groovy
roi = new Region("chr22:16591390-16591800")
println("Does variant $v overlap $roi? : " + roi.overlaps(v))
```
right:
**Output**
`````text
Does variant chr22:16591593 A/G overlap chr22:16591390-16591800? : true
`````
````````

# Querying VCFs for Presence of Variants

Finding out whether a particular variant is contained in a VCF is easy because the `VCF` class implements
the Groovy [membership operator](https://groovy-lang.org/operators.html#_membership_operator) which enables
use of the `in` keyword. Hence we can query if a VCF contains a particular variant by using code such as:

```groovy
vcf = VCF.parse('test.vcf')
variant = vcf[0]
if(variant in vcf)
    println "Variant $v is in the VCF!"
```

Note: in testing whether a variant is "in" a VCF, the position and allele are compared. However, the dosage (zygosity)
are not compared. Therefore, querying if a homozygous variant is present in a VCF that contains the VCF in heterozygous form
will return `true`. In fact, even if the variant is present but genotyped as homozygous reference, the query will return `true`.
To ascertain whether any sample has non-zero dosage for the variant, use the [find](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/VCF.html#find(gngs.Variant))
to retrieve the variant and then directly query the dosages. For example: `vcf.find(v).dosages.any { it > 0 }`

Warning: this function may behave in unexpected ways when multi-allelic variants are used as query or target VCF. The
`find`  and `in` functions will match a variant in a VCF if *any* allele matches *any* allele in the query.

# Info Fields

Info fields for each [Variant](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/Variant.html) are parsed and accessible as
a `Map` keyed on the name. Values are represented as `String`s and are *not* further parsed into typed entities matching the `INFO` fields declared in the
VCF header.

# Annotated Variants

The [Variant](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/Variant.html) class has built in support for parsing 
[VEP](https://asia.ensembl.org/info/docs/tools/vep/index.html) annotations. These are accessed via the 
`vepInfo` property. As VEP can create multiple annotations, the result is returned as a `List`, each entry which is
a `Map<String,Object>` of VEP annotations.


# Updating Variant Attributes

GNGS contains limited support for updating variants, allowing a variant that has been read from a VCF file to be modified
and then written out with the changes intact. Due to how `Variant` entries are stored, it is critical that when
updating a variant the update is performed within a `update` closure:


```groovy { title: "Updating a Variant to add an INFO property" }
v.update {
    v.info.MYINFO = "a new INFO property"
}
```

When updated this way, changes to INFO and genotype fields will be written out if made as part of a `filter` operation
or the VCF is written using the [print](https://ssadedin.github.io/groovy-ngs-utils/doc/gngs/VCF.html#print(java.lang.Appendable)) method.