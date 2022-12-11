// vim: ts=4:sw=4:expandtab:cindent:
/*
 * Copyright (c) 2012 MCRI, authors
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
package gngs.sample

import java.text.ParseException

import htsjdk.variant.vcf.VCFFileReader
import htsjdk.variant.vcf.VCFHeader;

import htsjdk.samtools.SamReader
import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.ValidationStringency

import com.xlson.groovycsv.CsvIterator;
import com.xlson.groovycsv.CsvParser;
import com.xlson.groovycsv.PropertyMapper


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
     *  - cram
     */
    Map    files = new Hashtable() // thread safe
	
	static List<String> SIMPLE_COLUMNS = [
        "Sample_ID",
        "Batch",
        "Cohort",
        "Fastq_Files",
        "Prioritised_Genes",
        "Sex",
        "Sample_Type",
        "Consanguinity",
        "Variants_File",
        "Pedigree_File",
        "Ethnicity",
        "VariantCall_Group",
        "DNA_Concentration",
        "DNA_Volume",
        "DNA_Quantity",
        "DNA_Quality",
        "DNA_Date",
        "Capture_Date",
        "Sequencing_Date",
        "Mean_Coverage",
        "Duplicate_Percentage",
        "Machine_ID",
        "Hospital_Centre",
        "Sequencing_Contact",
        "Pipeline_Contact",
        "Notes"]
    
    /**
     * MGHA redefined column order and contents to have a lot of things not of interest to others,
     * so have a separate mapping for them.
     */
    static List<String> MG_COLUMNS = [
		"Batch","Sample_ID","DNA_ID","Sex","DNA_Concentration","DNA_Volume","DNA_Quantity","DNA_Quality","DNA_Date","Cohort","Sample_Type",
        "Fastq_Files","Prioritised_Genes","Consanguinity","Variants_File",
        "Pedigree_File",
        "Ethnicity","VariantCall_Group","Capture_Date","Sequencing_Date","Mean_Coverage","Duplicate_Percentage","Machine_ID",
        "DNA_Extraction_Lab","Sequencing_Lab","Library_Preparation","Barcode_Pool_Size","Read_Type","Machine_Type","Sequencing_Chemistry",
        "Sequencing_Software","Demultiplex_Software","Hospital_Centre",
        "Hospital_Centre","Sequencing_Contact","Pipeline_Contact","Notes"
    ]
	
	/** Id of batch in which the sample was sequenced */
	String batch
	
	/** Whether the sample type is normal or tumor */
	SampleType sampleType = SampleType.NORMAL
	
	/** Whether the sample is consanguinous */
	Consanguinity consanguinity = Consanguinity.NOT_CONSANGUINEOUS

    /** Target (flagship) name */
    String  target

    /** List of genes prioritised for the sample */
    Map<String,Integer>    geneCategories

    /** The library */
    String library

    /** The sex of the sample */
    Sex sex = Sex.UNKNOWN
	
	Ethnicity ethnicity
	
    /** The pedigree of the family */
    String pedigree
	
	/** DNA quality in nanograms */
	float dnaConcentrationNg
	
	float dnaQuality
	
	float dnaQuantity
	
	List<Date> dnaDates
	
	List<Date> captureDates
	
	List<Date> sequencingDates
	
	/** Mean coverage as reported by sequencing provider */
	float meanCoverage
	
	/** Hospital or organization responsible for the patient from which the sample originated */
	String institution
	
	List<String> machineIds
	
	String sequencingContact
	
	String analysisContact
    
    String variantsFile
    
    /**
     * Recognised alternative identifiers for the sample
     */
    List<String> altIds

    Map<String,List<String>> fileMappings = [
        cram : ["cram"],
        bam : ["bam"],
        fastq: ["fastq.gz", "fq.gz"],
        coverage: ["exoncoverage.txt"],
        vcf : ["vcf"]
    ]

    void indexFileTypes() {
        fileMappings.each { index, endings -> indexFileType(index,endings) }
    }

    void indexFileType(String index, List<String> endings) {
        def matchingFiles  = files.all.grep { endings.any { String ending -> it.endsWith(ending) } }
        if(matchingFiles)
            this.files[index] = matchingFiles
    }

    /**
     * Parse the given file to extract sample info, where the file is in the
     * extended Melbourne Genomics Health Alliance format.
     * 
     * @param fileName
     * @return
     */
    static Map<String,SampleInfo> parse_mg_sample_info(String fileName) {
        SampleInfo.parse_sample_info(fileName, MG_COLUMNS)
    }
    
    /**
     * Parse sample info using auto-detection to determine if the format is
     * MGHA extended format or simplified format.
     * @param fileName
     * @return
     */
    static Map<String,SampleInfo> parse_sample_info(String fileName) {
        
        // We just read the first line and sniff the values
        List lines = readSampleInfoLines(fileName)
        if(lines.isEmpty())
            throw new RuntimeException("Sample information file is empty: please ensure at least one sample is specified in the file.")
           
        String [] line0 = lines[0].split('\t')
        
        int columns = line0.size()
        
        if(columns == MG_COLUMNS.size())
            return parse_mg_sample_info(fileName, MG_COLUMNS)
        else
        if(columns >= 5 && columns <= SIMPLE_COLUMNS.size())
            return parse_sample_info(fileName, SIMPLE_COLUMNS)
        else
            throw new RuntimeException(
                "Sample information file does not have expected number of columns. Expected either 5 - " +
                SIMPLE_COLUMNS.size() + " columns or exactly " + MG_COLUMNS.size() + " columns but observed $columns") 
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
    static Map<String,SampleInfo> parse_sample_info(String fileName, List columns) {
        
        String col0 = columns[0].toLowerCase()

        def lines = readSampleInfoLines(fileName)

        // Pad with optional blank fields
        lines = lines.collect { line ->
            def fields = line.split("\t")
            if(fields.size() > columns.size()) {
                throw new RuntimeException("Unable to parse sample information: file has more columns (${fields.size()}) than expected ($columns.size())")
            }
            return (fields + [""] * (columns.size() - fields.size())).join("\t")
        }

        int lineCount = 0
        CsvIterator parser = CsvParser.parseCsv(new StringReader(lines.join("\n")), columnNames: columns, readFirstLine: true, separator: '\t')
        def sample_info = parser.collect { PropertyMapper fields ->
            //				println "Found sample " + fields.Sample_ID

            String sampleId = fields.Sample_ID.tokenize(':')[0]
            
            List<String> altIds = fields.Identifiers.tokenize(':')

            try {
                def si = new SampleInfo(
                        sample: sampleId,
                        target: fields.Cohort,
                        sampleType: SampleType.decode(fields.Sample_Type),
                        geneCategories:  [:],
                        batch: fields.Batch,
                        pedigree: fields.Pedigree_File,
                        altIds: altIds
                        )

                si.sex = Sex.decode(fields.Sex)
                si.consanguinity = Consanguinity.decode(fields.Consanguinity)
                
                def if_field = { name, then_do ->
                    if(columns.contains(name) && fields[name])
                        then_do()
                }
                
                if_field('Ethnicity') {
                    si.ethnicity  = Ethnicity.decode(fields.Ethnicity)
                }
                
                if_field('DNA_Date') {
                    si.dnaDates = fields.DNA_Date?.split(",")*.trim().collect { parseDate(it) }
                }

                if_field('Capture_Date') {
                    si.captureDates = fields.Capture_Date.split(",")*.trim().collect { parseDate(it) }
                }

                if_field('Sequencing_Date') {
                    si.sequencingDates = fields.Sequencing_Date?.split(",")*.trim().collect { parseDate(it) }
                }

                if_field('DNA_Concentration') {
                    si.dnaConcentrationNg = fields.DNA_Concentration?.toFloat()
                }
                if_field('DNA_Quality') {
                    si.dnaQuality = fields.DNA_Quality?.toFloat()
                }
                if_field('DNA_Quantity') {
                    si.dnaQuantity = fields.DNA_Quantity?.toFloat()
                }
                if_field('Mean_Coverage') {
                    si.meanCoverage = fields.Mean_Coverage?.toFloat()
                }
                if_field('Variants_File') {
                    si.variantsFile = fields.Variants_File
                }
                if_field('Machine_ID') {
                    si.machineIds = fields.Machine_ID?.split(",")*.trim() as List
                }

                if_field('Sequencing_Contact') {
                    si.sequencingContact = fields.Sequencing_Contact
                }
                if_field('Pipeline_Contact') {
                    si.analysisContact = fields.Pipeline_Contact
                }
                if_field('Hospital_Centre') {
                    si.institution = fields.Hospital_Centre
                }
                si.files.all = fields.Fastq_Files.split(",")*.trim().collect {new File(it).parentFile?it:"../data/$it"}
                si.indexFileTypes()

                // Index category to gene
                if_field('Prioritised_Genes') {
                    def genes = fields.Prioritised_Genes.split(" ")*.split(":").collect {
                        // HACK: Excel is exporting weird white space characters that are not
                        // trimmed
                        int category = it[0].replaceAll("[^0-9]","").trim().toInteger()
                        [ category, /* genes */ it[1].split(",")*.trim() ]
                    }.collectEntries()

                    // Invert
                    genes.each { k,v -> v.each { si.geneCategories[it] = k } }
                }

                ++lineCount
                return si
            }
            catch(Exception e) {
                throw new RuntimeException("Error parsing meta data for sample ${fields.Sample_ID} on line $lineCount", e)
            }
        }.collectEntries { [it.sample, it] } // Convert to map keyed by sample name
    }

    static List<String> readSampleInfoLines(String fileName) {
        String simpleCol0 = SIMPLE_COLUMNS[0].toLowerCase()+"\t"
        String mgCol0 = MG_COLUMNS[0].toLowerCase() + "\t"
        new File(fileName).readLines().grep {
            !it.trim().startsWith('#') && // ignore comment lines
                    !it.trim().toLowerCase().startsWith(simpleCol0) && // ignore header line, if it is present
                    !it.trim().toLowerCase().startsWith(mgCol0) && // ignore header line, if it is present
                    it.trim() // ignore completely blank lines
        }
    }

    /**
     * Validates that text fields are in the correct format.
     * <p>
     * This format validation is specific to melbourne genomics.
     */
    void validate() {
        if(!(sample ==~ "[0-9]{9}"))
            throw new IllegalStateException("Sample ID ${sample} does not match prescribed format")

        if(!(batch ==~ "[0-9]{3}"))
            throw new IllegalStateException("Sample ID ${batch} does not match prescribed format")
    }

    /**
     * Return a tab separated string compatible with the samples.txt file format, 
     * containing the details for this sample.
     * <p>
     * 
     * @param father    optional father id, if provided, will be output into the pedigree field
     * @param mother    optional mother id, if provided, will be output into the pedigree field
     * 
     * @return
     */
    String toTsv(List<String> columns = SIMPLE_COLUMNS, String father=null, String mother=null) {
        /*
        [
            sample,
            batch,
            target,
            files.collect { it.key == "all" ? [] : it.value }.flatten().join(","),
            geneCategories.collect { it.key + ":" + it.value.join(",") }.join(" "),
            sex.encode()
        ].join("\t")
        */
        def mappings = [
       		"Batch" : { batch },
            "Sample_ID" : { sample },
            "Sex": { sex.encode() },
            "Cohort" : { target },
            "Fastq_Files": { files.collect { it.key == "all" ? [] : it.value }.flatten().join(",")},
            "Pedigree_File" : { 
                if(mother != null && father != null)  {
                    "familyID=$father,$mother"
                }
                else {
                    ""
                }
             }
        ]
        
        return columns.collect { key ->
            key in mappings ? mappings[key]() : ""
        }.join("\t")
    }

    String toString() {
        "$sample($sex)"
    }

    /**
     * Create a list of SampleInfo objects from provided files that
     * contain sample information in the header data. 
     * 
     * @param files
     * @return
     */
    static Map<String,SampleInfo> fromFiles(List<String> files, String mask=null) {
        // Convert to absolute path
        files = files.collect { new File(it).absolutePath }
        def collectBySample = { ext, extractSample ->

            def sampleExtracter = extractSample
            if(mask)
                sampleExtracter =  { f -> extractSample(f).replaceAll(mask,'') }

            def fs = files.grep {
                it.endsWith(ext)
            }

            return fs.groupBy(sampleExtracter).collectEntries {
                def f = [:]
                if(ext=="fastq.gz")
                    f["fastq"] = it.value
                else
                if(ext=="exoncoverage.txt" || ext=="coverage.txt" || ext=="cov.gz" || ext=="exoncoverage.txt.gz")
                    f["coverage"] = it.value
                else
                    f[ext] = it.value
                [it.key, [files: f, sample:it.key]]
            }
        }

        def samplesByFastq = collectBySample("fastq.gz") {
            new IlluminaFileNameParser(dialect:IlluminaFileNameParser.DIALECT.MGHA).parse(it).sample
        }

        def samplesByBam = collectBySample("bam") {
            getBAMSamples(it)[0].replaceAll('_$','') // legacy data had bad trailing _
        }

        def samplesByCram = collectBySample("cram") {
            getBAMSamples(it)[0].replaceAll('_$','') // legacy data had bad trailing _
        }

        def samplesByVcf = collectBySample("vcf") {
            getVCFSamples(it)[0].replaceAll('_$','') // legacy data had bad trailing _
        }

        def samplesByCoverage = collectBySample("exoncoverage.txt") {
            def sn = new File(it).name.replaceAll('\\..*$','')
            return sn
        }

        def samplesByGzCoverage = collectBySample("cov.gz") {
            def sn = new File(it).name.replaceAll('\\..*$','')
            return sn
        }

        def samplesByGzCoverage2 = collectBySample("exoncoverage.txt.gz") {
            def sn = new IlluminaFileNameParser().parse(it).sample
            return sn
        }

        def allFiles = [samplesByBam, samplesByCram, samplesByVcf, samplesByFastq, samplesByCoverage, samplesByGzCoverage, samplesByGzCoverage2]

        // Merge files from all of them
        def allSamples = allFiles*.keySet().sum().unique()
        Map<String,SampleInfo> result = allSamples.collect { s ->
            def sampleFiles = [all:[]]
            allFiles.each { samplesByType ->
                if(samplesByType[s]) {
                    sampleFiles += samplesByType[s].files
                    samplesByType[s].files.each { fileType, fileList ->
                        sampleFiles.all += fileList
                    }
                }
            }
            new SampleInfo(sample:s, files: sampleFiles)
        }.collectEntries { [ it.sample, it ] }

        return result
    }
    
    static List<String> getBAMSamples(String fileName) {
        SamReaderFactory samReaderFactory =
                SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)

        SamReader bam = samReaderFactory.open(new File(fileName))
        try {
            bam.getFileHeader().getReadGroups()*.sample
        }
        finally {
            bam.close()
        }
    }
    
    static List<String> getVCFSamples(String fileName) {
        VCFFileReader vcf = new VCFFileReader(new File(fileName), false)
        try {
            VCFHeader header = vcf.getFileHeader()
            return header.getSampleNamesInOrder()
        }
        finally {
            vcf.close()
        }
    }
        
    

    static Date parseDate(String dateValue) {
        if(!dateValue)
            return null
        return Date.parse("yyyyMMdd", dateValue)
    }
}
