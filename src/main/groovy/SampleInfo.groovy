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
import java.text.ParseException;

import com.xlson.groovycsv.CsvIterator;
import com.xlson.groovycsv.CsvParser;

enum Sex {
	MALE, FEMALE, OTHER, UNKNOWN
	
    private static codes = [
            1 : MALE,
            2 : FEMALE,
            "1" : MALE,
            "2" : FEMALE,
            "MALE" : MALE,
            "FEMALE": FEMALE,
            "Male" : MALE,
            "Female": FEMALE,
            "Unknown" : UNKNOWN,
            "other" : UNKNOWN // ped file compatibility
    ]
    
    String encode() {
        switch(this) {
            case MALE: "Male"; break
            case FEMALE: "Female"; break
            case OTHER: "Other"; break;
            case UNKNOWN: "Unknown"; break;
        }
    }
        
	static Sex decode(def  value) {
        if(value instanceof String)
            value = value?.trim()
		if(!value)
			return FEMALE
            
        if(codes.containsKey(value))
            return codes[value]

		throw new IllegalArgumentException("Bad sex value [$value] specified")
	}
}

enum SampleType {
	NORMAL, TUMOR
	
    private static Map codes = [
            "1" : NORMAL,
            "2" : TUMOR,
            "Normal" : NORMAL,
            "Tumour": TUMOR
    ]
 
	static SampleType decode(String value) {
        value = value?.trim()
		if(!value)
			return NORMAL
            
        if(codes.containsKey(value))
            return codes[value]
        
		throw new IllegalArgumentException("Bad sample type value [$value] specified")
	}
}

enum Consanguinity {
	NOT_CONSANGUINEOUS, CONSANGUINEOUS, SUSPECTED, UNKNOWN
    
    private static Map codes = [
            "0" : NOT_CONSANGUINEOUS,
            "1" : CONSANGUINEOUS,
            "2" : SUSPECTED,
            "8" : UNKNOWN,
            "No" : NOT_CONSANGUINEOUS,
            "Yes": CONSANGUINEOUS,
            "Suspected" : SUSPECTED,
            "Unknown" : UNKNOWN
    ]
	
	static Consanguinity decode(String value) {
        value = value?.trim()
        
        // Not strictly to spec: but this is the only non-optional field of many
        // so by itself it forces you to enter many other columns if it is required
        if(!value)
            return UNKNOWN
        
        if(codes.containsKey(value))
            return codes[value]
         
		throw new IllegalArgumentException("Bad consanguinity value [$value] specified")
	}
}

enum Ethnicity {
	UNKNOWN, EUROPEAN, AFRICAN, ASIAN
	
    private static Map codes = [
            "0" : UNKNOWN,
            "1" : EUROPEAN,
            "2" : AFRICAN,
            "European" : EUROPEAN,
            "Asian": ASIAN,
            "African" : AFRICAN,
            "Unknown" : UNKNOWN
    ]
    
	static Ethnicity decode(String value) {
        value = value?.trim()
		if(!value)
			return UNKNOWN
            
        if(codes.containsKey(value))
            return codes[value]
        
		throw new IllegalArgumentException("Bad ethnicity value [$value] specified")
	}
}

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
	
	static List<String> SIMPLE_COLUMNS = ["Sample_ID","Batch","Cohort","Fastq_Files","Prioritised_Genes","Sex","Sample_Type","Consanguinity","Variants_File","Pedigree_File","Ethnicity","VariantCall_Group","DNA_Concentration","DNA_Volume","DNA_Quantity","DNA_Quality","DNA_Date","Capture_Date","Sequencing_Date","Mean_Coverage","Duplicate_Percentage","Machine_ID","Hospital_Centre","Sequencing_Contact","Pipeline_Contact","Notes"]
    
    /**
     * MGHA redefined column order and contents to have a lot of things not of interest to others,
     * so have a separate mapping for them.
     */
    static List<String> MG_COLUMNS = [
		"Batch","Sample_ID","DNA_ID","Sex","DNA_Concentration","DNA_Volume","DNA_Quantity","DNA_Quality","DNA_Date","Cohort","Sample_Type",
        "Fastq_Files","Prioritised_Genes","Consanguinity","Variants_File",
        "Pedigree_File","Ethnicity","VariantCall_Group","Capture_Date","Sequencing_Date","Mean_Coverage","Duplicate_Percentage","Machine_ID",
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

    Map<String,String> fileMappings = [
        cram : "cram",
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
     * Parse the given file to extract sample info, where the file is in the
     * extended Melbourne Genomics Health Alliance format.
     * 
     * @param fileName
     * @return
     */
    static parse_mg_sample_info(fileName) {
        parse_sample_info(fileName, MG_COLUMNS)
    }
    
    /**
     * Parse sample info using auto-detection to determine if the format is
     * MGHA extended format or simplified format.
     * @param fileName
     * @return
     */
    static parse_sample_info(String fileName) {
        // We just read the first line and sniff the values
        List lines = readSampleInfoLines(fileName)
        if(lines.isEmpty())
            throw new RuntimeException("Sample information file is empty: please ensure at least one sample is specified in the file.")
           
        TSV iter = new TSV(new StringReader(lines[0]), readFirstLine:true)
        def row = iter[0]
        int columns = row.values.size()
        
        if(columns == MG_COLUMNS.size())
            return parse_sample_info(fileName, MG_COLUMNS)
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
    static parse_sample_info(fileName, columns) {
        
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
        def sample_info = parser.collect { fields ->
//				println "Found sample " + fields.Sample_ID
			
				try {
	                def si = new SampleInfo(
	                    sample: fields.Sample_ID,
	                    target: fields.Cohort,
						sampleType: SampleType.decode(fields.Sample_Type),
	                    geneCategories:  [:],
						batch: fields.Batch,
	                    pedigree: fields.Pedigree_File               
					)
					
					si.sex = Sex.decode(fields.Sex) 
					si.consanguinity = Consanguinity.decode(fields.Consanguinity)
					si.ethnicity  = Ethnicity.decode(fields.Ethnicity)
					if(fields.DNA_Date)
						si.dnaDates = fields.DNA_Date?.split(",")*.trim().collect { parseDate(it) }
					if(fields.Capture_Date)
						si.captureDates = fields.Capture_Date.split(",")*.trim().collect { parseDate(it) }
					if(fields.Sequencing_Date)
						si.sequencingDates = fields.Sequencing_Date?.split(",")*.trim().collect { parseDate(it) }
					
					if(fields.DNA_Concentration)
						si.dnaConcentrationNg = fields.DNA_Concentration?.toFloat()
					if(fields.DNA_Quality)
						si.dnaQuality = fields.DNA_Quality?.toFloat()
					if(fields.DNA_Quantity)
						si.dnaQuantity = fields.DNA_Quantity?.toFloat()
					if(fields.Mean_Coverage)
						si.meanCoverage = fields.Mean_Coverage?.toFloat()
					if(fields.Variants_File)
						si.variantsFile = fields.Variants_File
					if(fields.Machine_ID)
						si.machineIds = fields.Machine_ID?.split(",")*.trim() as List
						
					
					si.sequencingContact = fields.Sequencing_Contact
					si.analysisContact = fields.Pipeline_Contact
					si.institution = fields.Hospital_Centre
					
	                si.files.all = fields.Fastq_Files.split(",")*.trim().collect {new File(it).parentFile?it:"../data/$it"}
	                si.indexFileTypes()
	
	                // Index category to gene
	                if(fields.Prioritised_Genes) {
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
        String simpleCol0 = SIMPLE_COLUMNS[0].toLowerCase()
        String mgCol0 = MG_COLUMNS[0].toLowerCase()
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
    
    String toTsv() {
        [   
            sample, 
            batch, 
            target, 
            files.collect { it.key == "all" ? [] : it.value }.flatten().join(","), 
            geneCategories.collect { it.key + ":" + it.value.join(",") }.join(" "), 
            sex.encode()
        ].join("\t")
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
        files = files.collect { new File(it).canonicalFile.absolutePath }
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
                if(ext=="exoncoverage.txt" || ext=="coverage.txt" || ext=="cov.gz")
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
                new SAM(it).samples[0].replaceAll('_$','') // legacy data had bad trailing _
        }
        
        def samplesByCram = collectBySample("cram") {
                new SAM(it).samples[0].replaceAll('_$','') // legacy data had bad trailing _
        }
         
        def samplesByVcf = collectBySample("vcf") { 
            new VCF(it).samples[0].replaceAll('_$','') // legacy data had bad trailing _
        }
        
        def samplesByCoverage = collectBySample("exoncoverage.txt") { 
            def sn = new File(it).name.replaceAll('\\..*$','')
            return sn
        }
        
        def samplesByGzCoverage = collectBySample("cov.gz") { 
            def sn = new File(it).name.replaceAll('\\..*$','')
            return sn
        }
        
        def allFiles = [samplesByBam, samplesByCram, samplesByVcf, samplesByFastq, samplesByCoverage, samplesByGzCoverage] 
          
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
	
	static Date parseDate(String dateValue) {
		if(!dateValue)
			return null
		return Date.parse("yyyyMMdd", dateValue)
	}
}
