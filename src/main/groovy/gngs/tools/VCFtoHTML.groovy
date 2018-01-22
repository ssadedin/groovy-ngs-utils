/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2017 Simon Sadedin, ssadedin<at>gmail.com
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
package gngs.tools

import org.codehaus.groovy.runtime.StackTraceUtils;

import gngs.BED
import gngs.Cli
import gngs.Pedigrees
import gngs.ProgressCounter
import gngs.Regions
import gngs.SAM
import gngs.Utils
import gngs.VCF
import gngs.VCFSummaryStats
import gngs.Variant
import au.com.bytecode.opencsv.CSVWriter
import groovy.util.logging.Log
import groovy.xml.StreamingMarkupBuilder

import java.util.regex.Pattern

@Log
class VCFtoHTML {
    private boolean hasSNPEFF
    private boolean hasVEP
    private Map consColumns
    
    private Map baseColumns = new LinkedHashMap()
    private List<String> filters
    private BED target
    private Pedigrees pedigrees

    String banner = """
    """
    
    final static Pattern LOWER_CASE_BASE_PATTERN = ~/[agct]/
    
    OptionAccessor opts
    
    /**
     * Count of variants processed at any given time
     */
    int count = 1
    
    /**
     * Line that processing is up to in the VCF file
     */
    int lineIndex  
    
    int lastLines = 0
    
    Variant lastVariant  
    
    VCFSummaryStats stats = new VCFSummaryStats()
    
    List<String> exportSamples
    
    Map<String,SAM> bams = [:]
    
    float MAF_THRESHOLD=0.05
        
    VCFtoHTML(OptionAccessor opts) {
        this.opts = opts
        filters = []
        if(opts.filters) {
            filters = opts.filters
        }
        
        if(opts.maxMaf)
            MAF_THRESHOLD=opts.maxMaf.toFloat()
    }
    
    def js = [
        "http://ajax.aspnetcdn.com/ajax/jQuery/jquery-2.1.0.min.js",
        "http://cdn.datatables.net/1.10.0/js/jquery.dataTables.min.js",
        "http://ajax.googleapis.com/ajax/libs/jqueryui/1.10.4/jquery-ui.min.js",
        "http://cdnjs.cloudflare.com/ajax/libs/jquery-layout/1.3.0-rc-30.79/jquery.layout.min.js",
    //    "http://igv.org/web/beta/igv-beta.js"
    ]
        
    def css = [
        "http://cdn.datatables.net/1.10.0/css/jquery.dataTables.css",
        "http://ajax.googleapis.com/ajax/libs/jqueryui/1.10.4/themes/smoothness/jquery-ui.css",
        "http://cdnjs.cloudflare.com/ajax/libs/jquery-layout/1.3.0-rc-30.79/layout-default.css",
        "http://fonts.googleapis.com/css?family=PT+Sans:400,700",
        "http://fonts.googleapis.com/css?family=Open+Sans'",
        "http://maxcdn.bootstrapcdn.com/font-awesome/4.2.0/css/font-awesome.min.css",
    //     "http://igv.org/web/beta/igv-beta.css"
    ]
        
    def EXCLUDE_VEP = ["synonymous_variant","intron_variant","intergenic_variant","upstream_gene_variant","downstream_gene_variant","5_prime_UTR_variant"]
        
    def VEP_PRIORITY = [
            "transcript_ablation",
            "splice_donor_variant",
            "splice_acceptor_variant",
            "stop_gained",
            "frameshift_variant",
            "stop_lost",
            "initiator_codon_variant",
            "inframe_insertion",
            "inframe_deletion",
            "missense_variant",
            "transcript_amplification",
            "splice_region_variant",
            "incomplete_terminal_codon_variant",
            "synonymous_variant",
            "stop_retained_variant",
            "coding_sequence_variant",
            "mature_miRNA_variant",
            "5_prime_UTR_variant",
            "3_prime_UTR_variant",
            "non_coding_exon_variant",
            "nc_transcript_variant",
            "intron_variant",
            "NMD_transcript_variant",
            "upstream_gene_variant",
            "downstream_gene_variant",
            "TFBS_ablation",
            "TFBS_amplification",
            "TF_binding_site_variant",
            "regulatory_region_variant",
            "regulatory_region_ablation",
            "regulatory_region_amplification",
            "feature_elongation",
            "feature_truncation",
            "intergenic_variant"
    ]
    

    static void main(String [] args) {
        
        Utils.configureSimpleLogging()

        Cli cli = new Cli(usage:"VCFtoHTML <options>")
        cli.with {
            i 'vcf File', args: Cli.UNLIMITED, required:true
            p 'PED file describing relationships between the samples in the data', args:1
            o 'Output HTML file', args:1, required:true
            f 'Comma separated list of families to export', args:1
            a 'comma separated aliases for samples in input VCF at corresponding -i position', args: Cli.UNLIMITED
            target 'Exclude variants outside this region', args:1
            chr 'Confine analysis to chromosome', args:Cli.UNLIMITED
            diff 'Only output variants that are different between the samples'
            maxMaf 'Filter out variants above this MAF', args:1
            maxVep 'Only write a single line for each variant, containing the most severe consequence'
            tsv 'Output TSV format with the same variants', args:1
            filter 'A groovy expression to be used as a filter on variants before inclusion', args: Cli.UNLIMITED
            prefilter 'A groovy expression to be used as a pre-filter before any other steps', args: Cli.UNLIMITED
            nocomplex 'Ignore complex variants (regions where multiple variants overlap) in diff mode'
            nomasked 'Ignore variants in regions of the genome that are masked due to repeats'
            stats 'Write statistics to file', args:1
            mask 'Alias sample ids to first group of given regexp', args:1
            bam 'BAM File to compute coverage for variants if not in VCF', args:Cli.UNLIMITED, required:false
            roc 'Save a table of sensitivity and precision with sample <arg> as reference', args: 1, required:false
        }
        
        OptionAccessor opts = cli.parse(args)
        if(!opts)
            System.exit(1)
            
        new VCFtoHTML(opts).run()
    }
    
    void run() {
        
        List<String> allSamples = opts.is.collect { new VCF(it).samples }.flatten()
        
        if(opts.p) {
            pedigrees = Pedigrees.parse(opts.p)
        }
        else {
            // Find all the samples
            log.info "Found samples: " + allSamples
            pedigrees = Pedigrees.fromSingletons(allSamples)
        }
        
        if(opts.target) {
            target = new BED(opts.target).load()
        }
        
        resolveExportSamples(allSamples)
        
        List<VCF> vcfs = loadVCFs()
           
        List<Regions> vcfRegions
        if(opts.nocomplex)
            vcfRegions = vcfs*.toRegions()
            
        if(opts.bams)
            bams = opts.bams.collect { new SAM(it) }.collectEntries { [it.samples[0], it ] }
        
        // Adjust sample ids in case of duplicates (it is common to want to compare two VCFs for the same
        // sample, generated by different means)
        renameSampleIds(vcfs, exportSamples, pedigrees)
        
        log.info "Pedigree subjects are " + pedigrees.subjects.keySet()
        
        
        List aliases = []
        if(opts.as)
            aliases = opts['as'].collect { it.split "," }
        else
        if(opts.mask) {
            aliases = vcfs*.samples.flatten().collect { s ->
                def matches = (s =~ opts.mask)
                if(matches) {
                    // if starts with number, prefix with 'S' (required for javascript filtering)
                    def newSampleId = matches[0][1]
                    if(newSampleId ==~ /^[0-9].*/) {
                        'S' + newSampleId
                    }
                    else
                        newSampleId
                }
                else {
                    s // leave unchanged
                }
            }
        }
        
        // -------- Handle Aliasing of Samples ----------------
        if(aliases) {
            Map<String,String> sampleMap = [ vcfs*.samples.flatten(), aliases.flatten() ].transpose().collectEntries()
            
            log.info "Sample map = " + sampleMap
            
            vcfs.each { VCF vcf ->
                for(String s in vcf.samples) {
                    if(s in pedigrees.subjects.keySet()) { // may have been removed by filter on export samples
                        log.info "Pedigree rename $s => " + sampleMap[s]
                        pedigrees.renameSubject(s, sampleMap[s])
                    }
                    vcf.renameSample(s, sampleMap[s])
                }
            }
            
            exportSamples = exportSamples.collect { id -> sampleMap[id] }
        }
        
        log.info "Samples in vcfs are: " + vcfs.collect { vcf -> vcf.samples.join(",") }.join(" ")
        log.info "Export samples are: " + exportSamples
        
        checkAnnotations(vcfs)
        
        initColumnMappings()
        
        def tsvWriter = null
        if(opts.tsv)
            tsvWriter = new CSVWriter(new File(opts.tsv).newWriter(), '\t' as char)
        
            
        // Default to writing out the VCF file name replacing vcf extension with html
        String outputFileName = opts.o ?: opts.is[0].replaceAll('\\.vcf$','.html')
            
        new File(outputFileName).withWriter { w ->
           
            // Merge all the VCFs together
            VCF merged = vcfs[0]
            
            if(vcfs.size()>1) {
                vcfs[1..-1].eachWithIndex { vcf, vcfIndex ->
                    Utils.time("Merge VCF $vcfIndex") {
                         merged = merged.merge(vcf)
                    }
                }
            }
                
            log.info "Merged samples are " + merged.samples
            
            if(opts.tsv)
                tsvWriter.writeNext((baseColumns*.key+ consColumns*.key +exportSamples) as String[])
                    
            printHeadContent(w)
            
            lineIndex=0;
            
            lastVariant = null;
            Variant current  = null
            ProgressCounter progress = new ProgressCounter(withRate: true, withTime:true, extra: {
                stats.toString()
            })
            merged.each { v->
                try {
                    current = v
                    processVariant(v, w)
                }
                catch(Exception e) {
                    Exception e2 = new Exception("Failed to process variant " + current, e)
                    e2 = StackTraceUtils.sanitize(e2)
                    throw e2
                }
                progress.count()
            }
            
            progress.end()
            w.println """];"""
            
            w.println "var columnNames = ${json(baseColumns*.key + consColumns*.key + exportSamples)};"
            
            w.println """
            var samples = ${json(exportSamples)};
        
            var pedigrees = ${pedigrees.toJson()};
        
            var variantTable = null;
            \$(document).ready(function() {
                variantTable = \$.VariantTable('variantTable', samples, variants);
            });
            </script>
            <style type='text/css'>
            h1 {
                font-size: 16px;
            }
            #filterHelp {
                font-size: 9px;
            }
            label, .dataTables_info, .paginate_button, table#variantTable {
                font-size: 10px;
                font-family: verdana;
            }
        
            td.vcfcol { text-align: center; }
            tr.highlight, tr.highlighted, tr.highlight td.sorting_1 {
                background-color: #ffeeee !important;
            }
        
            #tags { float: right; }
                    
            .assignedTag, .rowTagDiv {
                background-color: red;
                border-radius: 5px;
                color: white;
                display: inline;
                padding: 3px 5px;
                font-size: 75%;
                margin: 3px;
            }
            .rowTagDiv { font-size: 55%; }
            .tag0 { background-color: #33aa33; }
            .tag1 { background-color: #3333ff; }
            .tag2 { background-color: #aa33aa; }
            .tag3 { background-color: #ff6600; }
            .unassignedTag { display: none; }
        
            </style>
            """;
            
           
            w.println """
            </head>
            <body>
            <p>Loading ...</p>
                <div class="ui-layout-north">
                    <h1>Variants ${opts.diff ? 'different between' : 'found in' } ${vcfs*.samples.flatten().unique().join(",")}</h1>
                    <span id=filterOuter><span id=filters> </span></span>
                    <div id=filterHelp style='display:none'> </div>
                </div>
                <div class="ui-layout-center">
                    <div id=tableHolder>
                    </div>
                </div>
                <div class="ui-layout-south">Report created ${new Date()}</div>
            </body>
            </html>
            """
        }
        
        if(tsvWriter != null)
            tsvWriter.close()
        
        printStats(stats)
    }

    private checkAnnotations(List vcfs) {
        hasVEP = true
        def noVeps = vcfs.findIndexValues { !it.hasInfo("CSQ") && !it.hasInfo("ANN") }
        if(noVeps) {
            System.err.println "INFO: This program requires that VCFs have VEP annotations for complete output. Output results will not have annotations and filtering\nmay be ineffective for samples in following files:"
            System.err.println "\n" + noVeps.collect { opts.is[(int)it]}.join("\n") + "\n"
            hasVEP = false
        }

        hasSNPEFF = vcfs.any { it.hasInfo("EFF") }
    }

    private initColumnMappings() {
        baseColumns += [
            'tags' : {''}, // reserved for tags
            'chr' : {it.chr },
            'pos' : {it.pos },
            'ref': {it.ref },
            'alt': {it.alt },
            'qual': {it.qual },
            'depth': {
                def dp = it.info.DP
                if(dp == null) {
                    bams*.value*.coverage(it.chr, it.pos).min()
                }
                else {
                    dp
                }
            },
            'families' : { v ->
                def fcount = v.pedigrees.count {  ped ->
                    def result = ped.samples.any {
                        v.sampleDosage(it)
                    }
                    return result
                }
                return fcount;
            }
        ]

        consColumns = [
            'gene' : {it['SYMBOL']},
            'cons' : {vep -> vep['Consequence'].split('&').min { VEP_PRIORITY.indexOf(it) } },
            'maf'  : this.&findMaxMaf
        ]
    }
    
    void processVariant(Variant v, Writer w) {
        
        ++stats.total
                    
        if((target != null) && !(v in target)) {
            ++stats.excludeByTarget
            return
        }
                    
        List baseInfo = baseColumns.collect { baseColumns[it.key](v) }
        List dosages = exportSamples.collect { v.sampleDosage(it) }
        if(dosages.every { it==0}) {
            ++stats.excludeNotPresent
            if(!opts.f) {
                println "WARNING: variant at $v.chr:$v.pos $v.ref/$v.alt in VCF but not genotyped as present for any sample"
            }
            return
        }
                        
        if(opts.diff && dosages.clone().unique().size()==1)  {
            ++stats.excludeByDiff
            return
        }
                        
        def refCount = v.getAlleleDepths(0)
        def altCount = v.getAlleleDepths(1)
                    
//        println v.toString() + v.line.split("\t")[8..-1] + " ==> " + "$refCount/$altCount"
                    
                    
        if(lineIndex++>0 && lastLines > 0)
            w.println ","
                    
        List<Object> consequences = [
            [
                Consequence: 'Unknown'
            ]
        ]
                    
        if(hasVEP) {
            try {
                if(opts.maxVep) {
                    consequences = [v.maxVep]
                }
                else {
                    consequences = v.vepInfo
                }
            }
            catch(Exception e) {
                // Ignore
            }
        }
        else
        if(hasSNPEFF) {
            consequences = v.snpEffInfo
        }
                    
        consequences = consequences.grep { it != null }
                    
        lastLines = 0
        boolean excludedByCons = true
        boolean printed = false
        consequences.grep { !it.Consequence.split('&').every { EXCLUDE_VEP.contains(it) } }.each { vep ->
                        
            if(EXCLUDE_VEP.contains(consColumns['cons'](vep))) {
                return
            }
            excludedByCons = false
                        
            if(hasVEP) {
                if(findMaxMaf(vep)>MAF_THRESHOLD)  {
                    ++stats.excludeByMaf
                    return
                }
            }
                        
            if(opts.nocomplex) {
                if(vcfRegions.any { Regions regions -> regions.getOverlaps(v.chr, v.pos-1,v.pos+1).size() > 1}) {
                    ++stats.excludeComplex
                    return
                }
            }
                        
            if(opts.nomasked) {
                if(v.alt.find(LOWER_CASE_BASE_PATTERN) || v.ref.find(LOWER_CASE_BASE_PATTERN)) {
                    ++stats.excludeByMasked
                    return
                }
            }
                        
            if(!filters.every { Eval.x(v, it) }) {
                ++stats.excludeByFilter
                return
            }
                        
            if(lastLines>0)
                w.println ","
                            
            List vepInfo = consColumns.collect { name, func -> func(vep) }
            w.print(groovy.json.JsonOutput.toJson(baseInfo+vepInfo+dosages + [ [refCount,altCount].transpose() ]))
            printed = true
                        
            if(opts.tsv)
                tsvWriter.writeNext((baseInfo+vepInfo+dosages) as String[])
            
            ++lastLines
        }
                    
        if(printed) {
            ++stats.totalIncluded
        }
                    
        lastVariant = v;
        if(excludedByCons)
            ++stats.excludeByCons
        
    }

    private printHeadContent(BufferedWriter w) {
        w.println """
            <html>
                <head>
            """
        w.println css.collect{"<link rel='stylesheet' href='$it'/>"}.join("\n").stripIndent()
        w.println js.collect{"<script type='text/javascript' src='$it'></script>"}.join("\n").stripIndent()

        // Embed the main vcf.js
        def vcfjs
        def fileVcfJsPath = "src/main/resources/vcf.js"
        if(new File(fileVcfJsPath).exists())
            vcfjs = new File(fileVcfJsPath).text
        else
            vcfjs = this.class.classLoader.getResourceAsStream("vcf.js").text

        w.println "<script type='text/javascript'>\n$vcfjs\n</script>"

        w.println "<script type='text/javascript'>"

        w.println "var variants = ["
    }

    void resolveExportSamples(List allSamples) {
        exportSamples = pedigrees.families.values().collect { it.individuals*.id }.flatten().grep { it in allSamples };

        def exportFamilies = pedigrees.families.keySet()
        if(opts.f) {
            exportFamilies = opts.f.split(",").collect { it.trim() }
            exportSamples = exportSamples.grep { s -> exportFamilies.any { f-> pedigrees.families[f].samples.contains(s)  }}
            pedigrees.families.keySet().grep { !exportFamilies.contains(it) }.each { pedigrees.removeFamily(it) }
        }
        log.info "Samples to export: $exportSamples"

        if(exportSamples.empty) {
            System.err.println "ERROR: No samples from families $exportFamilies found in VCF file $opts.i"
            System.err.println "\nSamples found in VCF are: " + (new VCF(opts.i)).samples
        }
    }

    List<VCF> loadVCFs() {
        
        List<String> preFilters = []
        if(opts.prefilters) {
            preFilters = opts.prefilters
        }
        
        List<VCF> vcfs = opts.is.collect { String vcfPath ->
            
            log.info "Read $vcfPath ..."
            
            VCF vcf = new VCF(vcfPath)
            List applicableExportSamples = exportSamples.grep { 
                it in vcf.samples
            }
            
            log.info "The following export samples are applicable to $vcfPath: ${applicableExportSamples.join(',')}"
            
            vcf = VCF.parse(new File(vcfPath), null, samples:applicableExportSamples?:null) { Variant  v ->
                
                if(opts.chrs && !(v.chr in opts.chrs ))
                    return false

                if(!preFilters.every { Eval.x(v, it) }) {
                    ++stats.excludeByPreFilter
                    return false
                }
            }
            
            log.info "Retained ${vcf.size()} variants from $vcfPath"
            return vcf
        }
        return vcfs
    }
    
    private void printStats(VCFSummaryStats stats) {
        println " Summary ".center(80,"=")
        
        Writer statsWriter = null
        if(opts.stats)
            statsWriter = new File(opts.stats).newWriter()
            
        try {    
            ["total", "excludeByPreFilter","excludeByDiff", "excludeByCons", "excludeByMaf", "excludeComplex", "excludeByMasked", "excludeByFilter", "excludeByTarget", "excludeNotPresent","totalIncluded"].each { prop ->
                println prop.padLeft(20) + " :" + stats[prop]
                if(opts.stats) {
                    statsWriter.println([prop,stats[prop]].join('\t'))
                }
            }
        }
        finally {
            if(statsWriter)
                statsWriter.close()
        }
    }

    /**
     *  If VCFs have samples with the same id, we now need to rename them to get a sensible result
     *  We also have to rename any instance of those in the "exportSamples" and the associated pedigrees
     *  since all the sample ids have to match for the downstream logic to work.
     * 
     * @param vcfs
     * @param exportSamples
     * @param pedigrees
     */
    private void renameSampleIds(List vcfs, List exportSamples, Pedigrees pedigrees) {
        List accumulatedSamples = vcfs[0].samples.clone()
        int n = 0
        if(vcfs.size() > 1) {
            for(VCF vcf in vcfs[1..-1]) {
                vcf.samples.collect { s ->
                    String newSample = s
                    int i = 1
                    while(newSample in accumulatedSamples) {
                        newSample = s + "_$i"
                        ++i
                    }
                    if((newSample != s) && (s in exportSamples)) {
                        println "Rename sample $s to $newSample in vcf $n"
                        exportSamples << newSample
                        if(pedigrees.subjects[s]) {
                            pedigrees.subjects[s].copySubject(s, newSample)
                            pedigrees.subjects[newSample] = pedigrees.subjects[s]
                        }

                        vcf.renameSample(s, newSample)
                    }
                    accumulatedSamples << newSample
                    return newSample
                }
                ++n
            }
        }
    }
    
    def json(obj) {
        groovy.json.JsonOutput.toJson(obj)
    }
    
    def findMaxMaf(vep) {
        [vep.EA_MAF, vep.ASN_MAF, vep.EUR_MAF, vep.ExAC_MAF].collect{ mafValue ->
            mafValue?mafValue.split('&'):[]
        }.flatten().collect {
            it.tokenize(':')
        }.grep {
            !it.isEmpty()
        }.collect {
            it[-1].toFloat()
        }.max() ?: 0.0f
    }
}
    
    
    
