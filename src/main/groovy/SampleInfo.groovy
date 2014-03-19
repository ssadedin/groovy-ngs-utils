// vim: ts=4:sw=4:expandtab:cindent:
/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2014 Simon Sadedin, ssadedin<at>gmail.com
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
/**
 * Meta data about a sample.
 * <p>
 * Designed to be compatible with the MGHA sample information format.
 */
class SampleInfo {

    /** Sample name */
    String  sample

    /** 
     * List of files containing data specific to this sample, indexed by file types:
     *  - fastq
     *  - coverage (output from coverageBed)
     *  - vcf
     *  - bam
     */
    Map    files = new Hashtable() // thread safe

    /** Target (flagship) name */
    String  target

    /** List of genes prioritised for the sample */
    Map<String,Integer>    geneCategories

    /** The library */
    String library

    /** The sex of the sample */
    String sex

    /** The pedigree of the family */
    String pedigree

    Map<String,String> fileMappings = [
        bam : "bam",
        fastq: "fastq.gz",
        coverage: "exoncoverage.txt",
        vcf : "vcf"
    ]

    void indexFileTypes() {
        fileMappings.each { index, ending -> indexFileType(index,ending) }
    }

    void indexFileType(String index, String ending) {
        def matchingFiles  = files.all.grep { it.endsWith(ending) }
        if(matchingFiles)
            this.files[index] = matchingFiles
    }

    /**
     * Parse the given file to extract sample information
     *
     * @return  List of (Map) objects defining properties of each
     *          sample: 
     *              <li>sample name (sample), 
     *              <li>FastQ files (files), 
     *              <li>Name of flagship (target)
     *              <li>Genes to be classed as high priority (genes)
     */
    static parse_sample_info(fileName) {
        def lines = new File(fileName).readLines().grep { !it.trim().startsWith('#') }
        def sample_info = lines.collect { it.split("\t") }.collect { fields ->
                fields = fields as List
                def si = new SampleInfo(
                    sample: fields[0], 
                    target: fields[1], 
                    geneCategories:  [:],
                    library: fields[0],
                    sex: fields[4],
                    pedigree: fields[8]
                ) 
                si.files.all = fields[2].split(",")*.trim().collect {new File(it).parentFile?it:"../data/$it"}
                si.indexFileTypes()

                if(fields.size()>3) {
                    fields[3].split(",")*.trim()

                    // Index category to gene
                    if(fields[3]) {
                        def genes = fields[3].split(" ")*.split(":").collect { 
                            [ /* category */ it[0].trim().toInteger(), /* genes */ it[1].split(",")*.trim() ]
                        }.collectEntries()

                        // Invert 
                        genes.each { k,v -> v.each { si.geneCategories[it] = k } }
                    }
                }

                return si
        }.collectEntries { [it.sample, it] } // Convert to map keyed by sample name
    }
    
    String toTsv() {
        [sample, target, files.collect { it.key == "all" ? [] : it.value }.flatten().join(","), geneCategories.collect { it.key + ":" + it.value.join(",") }.join(" "), sex].join("\t")
    }

    String toString() {
        "$sample($sex)"
    }
}
