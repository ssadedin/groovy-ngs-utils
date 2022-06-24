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

import gngs.VepConsequence
import org.codehaus.groovy.runtime.StackTraceUtils;

import gngs.*
import gngs.plot.Lines
import gngs.plot.Plot
import graxxia.Stats
import com.opencsv.CSVWriter
import de.erichseifert.gral.graphics.Drawable
import de.erichseifert.gral.ui.InteractivePanel
import groovy.json.JsonOutput
import groovy.transform.AutoClone
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import groovy.xml.StreamingMarkupBuilder

import java.text.NumberFormat
import java.util.regex.Pattern

import javax.swing.JFrame

class CustomInfo {
    String infoField
    String type
    String rendering
}

@CompileStatic
@AutoClone
class ROCPoint {
    
   
    Variant variant
    
    boolean truePositive
    
    int tpCount
    
    int fpCount
    
    int metricBin
}

abstract class MetricBins {
    String name 
    
    MetricBins(String name) {
        this.name = name
    }
    
    abstract int getBinIndex(Variant v, String sample) 
}

class DPMetricBins extends MetricBins {
    
    DPMetricBins() {
        super('Depth')
    }
    
    int getBinIndex(Variant v, String sample) {
        return v.getTotalDepth(sample)
    }
}

class QualMetricBins extends MetricBins {
    QualMetricBins () {
        super('Qual')
    }
    
    int getBinIndex(Variant v, String sample) {
        return Math.round(v.qual)
    }    
}

/**
 * A comprehensive tools supporting many different functions to render a VCF in HTML format,
 * as well as functionality to "diff" VCFs including computing ROC-like statistics.
 * 
 * @author Simon Sadedin
 */
@Log
class VCFtoHTML {
    private Map<String,String> sampleMap
    private String outputFileName
    private List<String> rocSamples
    private boolean hasSNPEFF
    private boolean hasVEP
    private Map consColumns
    
    private Map baseColumns = new LinkedHashMap()
    private List<String> filters
    
    private List<BED> targets = []
    
    private BED excludeRegions
    private Pedigrees pedigrees
    
    private List<CustomInfo> customInfos = []
    
    private Map<String,String> htmlFilters = [:]

    private NumberFormat THREE_DIGIT_PRECISION 
    
    
    {
        THREE_DIGIT_PRECISION = NumberFormat.getNumberInstance()
        THREE_DIGIT_PRECISION.maximumFractionDigits = 3
    }
    
    static final String DEFAULT_STYLES = """
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
            td.vcfcol a.igvLink.wasVisited { color: purple !important; }

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
            .rowTagDiv { font-size: 75%; }
            .tag0 { background-color: #33aa33; }
            .tag1 { background-color: #3333ff; }
            .tag2 { background-color: #aa33aa; }
            .tag3 { background-color: #ff6600; }
            .unassignedTag { display: none; }

            td.priority-1  { font-weight: bold; color: #ffaa00 !important; }
            td.priority-2  { font-weight: bold; color: #ffaa00 !important; }
            td.priority-3  { font-weight: bold; color: #ff6600 !important; }
            td.priority-4  { font-weight: bold; color: #ff0000 !important; }
            div.reportDesc { font-weight: normal; font-size: 92%; margin-bottom: 5px; }
        </style>
    """.stripIndent()

    String banner = """
    """
    
    final static Pattern LOWER_CASE_BASE_PATTERN = ~/[agct]/
    
    CliOptions opts
    
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
    
    float maxMaf=0.05
    
    List<ROCPoint> rocPoints = []
    
    Map<String, List<ROCPoint>> sampleROCPoints 
    
    /**
     * The number of variants included from the reference sample 
     * when roc option enabled
     */
    int rocTotalVariants = 0
    
    /**
     * The sample selected as the refernece sample when ROC mode is enabled
     */
    String referenceSample
    
    boolean allCons = false
    
    /**
     * If enabled, a writer that outputs differences as tab separated format
     */
    CSVWriter tsvWriter = null
    
    MetricBins bins = new DPMetricBins()
    
    double minROCSensitivity = 0.2
        
    FASTA fasta
    
    List<Closure> parsedFilters = []
    
    List<Closure> parsedPreFilters = []
    
    VCFtoHTML(CliOptions opts) {
        this.opts = opts
        filters = []
        if(opts.filters) {
            filters = opts.filters
            parsedFilters = filters.collect { 
                new GroovyShell().evaluate( "{ x ->\n\n$it\n}")
            }
        }
        
        if(opts.ref)
            this.fasta = new FASTA(opts.ref)
            
        if(opts.norep && !opts.ref)
            throw new IllegalArgumentException('If the norep argument is given, please also specify the reference using the "ref" option')
        
        if(opts.maxMaf)
            maxMaf=opts.maxMaf.toFloat()
            
        if(this.opts.rocBinSize)
            this.minBinSize = this.opts.rocBinSize.toInteger()
            
            
        if(this.opts.rocmetric) {
//                = new QualMetricBins()
            switch(this.opts.rocmetric) {
                case 'qual':
                    this.bins = new QualMetricBins()
                    break
                case 'depth':
                    this.bins = new DPMetricBins()
                    break
                default:
                    throw new IllegalArgumentException("Please use either 'depth' or 'qual' as the setting for ROC metric")
            }
        }
        
        if(this.opts.rocminsens) {
            this.minROCSensitivity = this.opts.rocminsens.toDouble()
        }
        
        if(this.opts.infos) {
            this.customInfos = this.opts.infos.collect {
                List parts = it.tokenize(':')
                
                CustomInfo info = new CustomInfo(infoField:parts[0])
                if(parts.size()>1)
                    info.type = parts[1]
                if(parts.size()>2)
                    info.rendering = parts[2]
                return info
            }
        }
        
        if(this.opts.htmlFilters) {
            this.htmlFilters = this.opts.htmlFilters.collectEntries { it.tokenize(':') }
        }
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
        
    List<String> excludedVEPConsequences = ["synonymous_variant","intron_variant","intergenic_variant","upstream_gene_variant","downstream_gene_variant","5_prime_UTR_variant"]
        
    static void main(String [] args) {
        
        Utils.configureSimpleLogging()

        Cli cli = new Cli(usage:"VCFtoHTML <options>")
        cli.with {
            i 'vcf File', args: Cli.UNLIMITED, required:true
            p 'PED file describing relationships between the samples in the data', args:1
            o 'Output HTML file', args:1, required:true
            f 'Comma separated list of families to export', args:1
            a 'comma separated aliases for samples in input VCF at corresponding -i position', args: Cli.UNLIMITED
            id 'Include variant ids in output'
            ref 'Reference FASTA to look up repeat regions, required for norep option', args: 1
            target 'Variants must fall within this region to be included (multiple allowed)', args:Cli.UNLIMITED
            xtarget 'Exclude variants inside these target regions from the comparison', args:1
            genelist 'Add gene priorities based on two column, tab separated file', args:1
            chr 'Confine analysis to chromosome', args:Cli.UNLIMITED
            diff 'Only output variants that are different between the samples'
            maxMaf 'Filter out variants above this MAF', args:1
            maxVep 'Only write a single line for each variant, containing the most severe consequence'
            allCons 'Do not filter by consequence (default: only variants with significant consequence included)'
            tsv 'Output TSV format with the same variants', args:1
            htmlFilter 'An expression to add automatically to the HTML report as a filter', args: '*'
            filter 'A groovy expression to be used as a filter on variants before inclusion', args: Cli.UNLIMITED
            prefilter 'A groovy expression to be used as a pre-filter before any other steps', args: Cli.UNLIMITED
            nocomplex 'Ignore complex variants (regions where multiple variants overlap) in diff mode'
            nomasked 'Ignore variants in regions of the genome that are masked due to repeats'
            stats 'Write statistics to file', args:1
            mask 'Alias sample ids to first group of given regexp', args:1
            bam 'BAM File to compute coverage for variants if not in VCF', args:Cli.UNLIMITED, required:false
            roc 'Save a table of sensitivity and precision with sample <arg> as reference', args: 1, required:false
            rocBinSize 'Minimum number of values to include in each ROC point', args:1, required: false
            rocCsv 'Save ROC data in CSV format in the given file', args:1, required:false
            rocmetric 'Metric to use for calculating ROC (qual, depth)', args:1, required: false
            rocminsens 'Minimum sensitivity to plot for ROC (20%)', args:1, required: false
            title 'Title to apply in various reports that are output', args:1, required: false
            desc 'Description to show in plain text under the title', args:1, required: false
            norep 'Do not include variants in or adjacent to repeat regions', required: false
            vcfjs 'Link to vcf javascript library instead of embedding', args:1, required: false, type: File
            info 'Custom INFO field to include in the report (provide multiple times)', args:'*', required: false
            hide 'Hide named columns from output in HTML', args: '*', required: false
            igvPadding 'Padding to apply to IGV links in bp', args: 1, required: false, type: Integer
        }
        
        CliOptions opts = new CliOptions(opts:cli.parse(args))
        if(!opts || !opts.i)
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
        
        if(opts.targets) {
            targets = opts.targets.collect { t ->
                log.info "Loading target region $t"
                new BED(t).load() 
            }
        }
        
        if(opts.xtarget) {
            excludeRegions = new BED(opts.xtarget).load()
        }
        
        resolveExportSamples(allSamples)
        
        if(opts.roc)
            referenceSample = opts.roc
        
        // Samples to perform ROC on are those not specified as the reference sample
        rocSamples = exportSamples.grep { it != referenceSample }
        
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
            aliases = applySampleMask(aliases, vcfs)
        }
        
        if(opts.allCons) {
            this.excludedVEPConsequences = []
        }
        
        // -------- Handle Aliasing of Samples ----------------
        if(aliases) {
            
            log.info "Samples from VCFs are: " + vcfs*.samples.flatten()
            
            sampleMap = [ vcfs*.samples.flatten(), aliases.flatten() ].transpose().collectEntries()
            
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
            rocSamples = exportSamples.grep { it != referenceSample }
        }
        else {
            sampleMap = vcfs*.samples.flatten().collectEntries { [it,it] }
        }
        
        if(opts.roc) {
            // Initialise with a list for each sample
            sampleROCPoints = rocSamples.collectEntries { [it, []] }
        }
        
        
        log.info "Samples in vcfs are: " + vcfs.collect { vcf -> vcf.samples?.join(",") }.join(" ")
        log.info "Export samples are: " + exportSamples
        
        checkAnnotations(vcfs)
        
        initColumnMappings()
        
        if(opts.tsv)
            tsvWriter = new CSVWriter(new File(opts.tsv).newWriter(), '\t' as char)
        
            
        // Default to writing out the VCF file name replacing vcf extension with html
        outputFileName = opts.o ?: opts.is[0].replaceAll('\\.vcf$','.html')
            
        new File(outputFileName).withWriter { w ->
           
            processVariantsAndWriteOutput(vcfs, vcfRegions, w)
        }
        
        if(tsvWriter != null)
            tsvWriter.close()
        
        printStats(stats)
        
        printROC()
    }

    private processVariantsAndWriteOutput(List vcfs, List vcfRegions, BufferedWriter w) {
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
                processVariant(v, vcfRegions, w)
            }
            catch(Exception e) {
                Exception e2 = new Exception("Failed to process variant " + current, e)
                e2 = StackTraceUtils.sanitize(e2)
                throw e2
            }
            progress.count()
        }

        progress.end()


        Map<String,Integer> genePriorities = [:]
        if(opts.genelist)
            genePriorities = new File(opts.genelist).readLines()*.tokenize('\t').collectEntries { [it[0],it[1].toInteger()] }

        w.println """];"""

        w.println "var columnNames = ${json(baseColumns*.key + consColumns*.key + customInfos*.infoField + exportSamples )};"
        w.println "var hiddenColumns = ${json(opts.hides?:['size'])};\n</script>"
        if(opts.title)
            w.println "<script type='text/javascript'>var customTitle = ${json(opts.title)};\n document.title=customTitle;</script>"

        w.println "<script type='text/javascript'>var customDesc = ${opts.desc?json(opts.desc):'null'};</script>"
        w.println "<script type='text/javascript'>var igvPadding = ${opts.igvPadding?opts.igvPadding:-1};</script>"

        injectVCFLibrary(w)

        w.println """
            <script type='text/javascript'>
            var samples = ${json(exportSamples)};

            var pedigrees = ${pedigrees.toJson()};

            var genePriorities = ${json(genePriorities)};
        
            var variantTable = null;
            \$(document).ready(function() {
                variantTable = \$.VariantTable('variantTable', samples, variants);
            });
            </script>

            $DEFAULT_STYLES
        """;

        
        def reportTitle =  opts.title ?: """Variants ${opts.diff ? 'different between' : 'found in' } ${vcfs*.samples.flatten().unique().join(",")}"""
        
        def descElement = opts.desc ? "<div class='reportDesc'>$opts.desc</div>" : ""

        w.println """
            </head>
            <body>
            <p>Loading ...</p>
                <div class="ui-layout-north">
                    <h1>${reportTitle}</h1>
                    ${descElement}
                    <span id=filterOuter><span id=filterButtons></span>&nbsp;<span id=filters> </span></span>
                    <div id=filterHelp style='display:none'> </div>
                </div>
                <div class="ui-layout-center">
                    <div id=tableHolder>
                    </div>
                </div>
                <div id='south-pane' class="ui-layout-south">
                    <div id='southContent' style='height: 100%'></div>
                    <div style='font-size:80%; margin-top: 4px;'>Report created ${new Date()}</div>
                </div>
                <iframe style='display:none' src='about:blank' id='igvframe' name='_igv'></iframe>
            </body>
            </html>
            """
    }

    private List applySampleMask(List aliases, List vcfs) {
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
        return aliases
    }
    
    /**
     * Calculate
     */
    void printROC() {
        
        if(rocTotalVariants == 0) {
            log.warning("No variants in reference sample: ROC cannot be computed")
            return
        }
        
        log.info "ROC table is based on sensitivity against ${rocTotalVariants} reference variants"
        
        Map<String,List<Map>> rocs = [:]

        for(rocSample in rocSamples) {
            
            List<Map> fullROC = calculateROC(rocSample)
            
            rocs[rocSample] = fullROC
            
            List<Map> roc = fullROC.grep { it.sensitivity > minROCSensitivity }
            
            File rocFile = new File(outputFileName.replaceAll('.html$', rocSample+'.roc.md'))
            rocFile.withWriter { w ->
                
                String hd = "ROC Table for ${rocSample} vs ${referenceSample}"
                w.println hd
                w.println ("=" * hd.size())
                w.println ""
                
                Utils.table(out:w,
                    [bins.name, 'Variants', 'Sensitivity', 'Precision'],
                    [roc*.bin, roc*.total, roc*.sensitivity.collect {Utils.perc(it)}, roc*.precision.collect{Utils.perc(it)}].transpose()
                )
            }
            
            if(opts.rocCsv) {
                saveROCCSV(bins, roc) 
            }
        }
        
        Plot rocPlot = new Plot(
            title: 'ROC Curve for ' + rocSamples.join(',') + ' (Ranked by ' + bins.name + ')', xLabel: 'False Positive Rate (1 - Precision)', 
            yLabel: 'Sensitivity') 
        
        
        rocs.each { rocSample, roc ->
            rocPlot << new Lines(x: roc*.precision.collect { (1 - it)*100 }, y: roc*.sensitivity*.multiply(100), displayName: rocSample) 
        }
            
        String rocImage = outputFileName.replaceAll('.html$', '.roc.png')
        rocPlot.save(rocImage)
        log.info "ROC image saved as ${rocImage}"
        
        String rocHTML = rocImage.replaceAll('\\.png', '\\.html')
        
        new File(rocHTML).text = getROCHTML(bins, rocs)
        
        // Save roc HTML
        log.info "ROC html saved as ${rocHTML}"
        
    }
    
    String getROCHTML(MetricBins bins, Map<String, List<Map>> points) {
        
        String titleExtra = opts.title ? "($opts.title)" :  ""
        String rocJson = JsonOutput.prettyPrint(JsonOutput.toJson(points.collect { [ name: it.key, roc: it.value ]}))
        return """<!DOCTYPE html>\n<html>
            <head>
                <script type='text/javascript'>
                var rocData = ${rocJson};
                </script>
                <script src='http://nvd3.org/assets/lib/d3.v3.js'></script>
                <script src='http://nvd3.org/assets/js/nv.d3.js'></script>
                <link href="http://nvd3.org/assets/css/nv.d3.css" rel="stylesheet">
            </head>
            <body>
                <h1>ROC Curve: Sensitivity and False Positive Rate by ${bins.name} $titleExtra</h1>
            
                <b>Max Metric Value:</b> <input type=text name=maxMetric id=maxMetric> <button onclick='redraw()'>Redraw</button>
            
                <div id=chart style='width: 800; height: 600;'>
                    <svg></svg>
                </div>
                <script type='text/javascript'>
            
            
                function redraw() {
                      document.getElementById('chart').innerHTML='<svg></svg>';
                      chart = nv.models.lineChart()
                                    .margin({left: 100})  
                                    .transitionDuration(350)  
                                    .showLegend(true)   
                                    .showYAxis(true)        
                                    .showXAxis(true)       
                                    // .interpolate('basis')
                      ;
            
                      chart.xAxis     //Chart x-axis settings
                          .axisLabel('False Positive Rate (1 - Precision)')
                          .tickFormat(d3.format(',r'));
            
                      chart.yAxis     //Chart y-axis settings
                          .axisLabel('Sensitvity')
                          .tickFormat(d3.format('.02f'));
            
            
                      let maxMetric = 999999999999999;
                      let maxMetricValue = document.getElementById('maxMetric').value
                      if(maxMetricValue)
                          maxMetric = parseFloat(maxMetricValue);

                      var myData  = rocData.map(rocObj => { return {
                         key: rocObj.name,
                         values: rocObj.roc.filter(point => point.bin < maxMetric)
                                       .map(point => { return { x: 1 - point.precision, y: point.sensitivity, metric: point.bin } })
                      }})
            
                      chart.forceX(0)

                      if(!maxMetricValue)
                          chart.forceY(0)
            
                      chart.tooltipContent((key, x, y, e, graph) => { 
                         let metric = e.point.metric
                          return '<b>Sensitivity:</b> ' + y + '<br><b>FPR:</b>' + x + '<br><b>Metric:</b>' + metric;
                      })
            
                      d3.select('#chart svg') 
                          .datum(myData)      
                          .call(chart);       
            
                      
                      console.log("Done")
                }
                redraw()
                </script>
            
            </body>
            </html>
        """.stripIndent()
    }
  
    
    void saveROCCSV(MetricBins bins, List<Map> points) {
        new File(opts.rocCsv).withWriter { w ->
            w.write([bins.name, 'Values', 'Sensitivity', 'Precision'].join(',') + '\n')
            for(rocPoint in points) {
                w.write([rocPoint.bin, rocPoint.total, rocPoint.sensitivity, rocPoint.precision].join(',') + '\n')
            }
        }
    }
    
    int minBinSize = 10
    
    /**
     * Compute an ROC style curve for the sample 
     * 
     * @param sample
     * @return
     */
    List<Map> calculateROC(String sample) {
        
        List<ROCPoint> rocPoints = sampleROCPoints[sample]
        if(!rocPoints) {
            log.info "No points found on the ROC curve for $sample: will exclude from output"
            return []
        }
       
        rocPoints.each { 
            it.metricBin = bins.getBinIndex(it.variant, sample)
        }
        
        // Sort from highest value of metric to lowest
        List<ROCPoint> sorted = 
            rocPoints.sort { -it.metricBin }
        
        sorted.inject { ROCPoint last, ROCPoint next ->
            next.tpCount = last.tpCount
            next.fpCount = last.fpCount
            if(next.truePositive) {
                ++next.tpCount 
            }
            else {
                ++next.fpCount
            }
            return next
        }
        
        TreeMap<Integer,List<ROCPoint>> rawRoc = sorted.groupBy { it.metricBin } 
        
        TreeMap<Integer,List<ROCPoint>> roc = new TreeMap<Integer, List<ROCPoint>>()
        
        List<ROCPoint> currentGroup = []
//        double binTotal = 0d
        rawRoc.each { Integer binIndex, List<ROCPoint> values ->
            currentGroup.addAll(values)
            if(currentGroup.size() > minBinSize) {
                // log.info "Merging " + currentGroup*.metricBin.unique() + " bins together (" + currentGroup*.metricBin.unique().join(',') + ") to achieve $minBinSize values " 
                roc.put(Math.round(Stats.mean(values*.metricBin)), currentGroup)
                currentGroup = []
            }
        }
        
        if(roc.isEmpty()) 
            throw new IllegalStateException("Too few points in each bin to create ROC data. Please use -rocBinSize to specify smaller bin size")
        
//        if(currentGroup.size() > minBinSize) {
//            roc.put(Math.round(Stats.mean(values*.metricBin)), currentGroup)
//        }
        
        return roc.collect { Map.Entry<Integer, List<ROCPoint>> binEntry ->
            
            List<ROCPoint> bin = binEntry.value
            
            double sensitivityEstimate = graxxia.Stats.mean(bin.collect { it.tpCount / rocTotalVariants })
            double precisionEstimate = 1.0d - graxxia.Stats.mean(bin.collect { it.fpCount / (1 + it.fpCount + it.tpCount) })
            
            [
                bin: binEntry.key,
                total: bin.size(), 
                binPrecision: bin.count { it.truePositive } / bin.size(),
                sensitivity: sensitivityEstimate,
                precision: precisionEstimate
            ]
        }
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
        ]
        
        if(opts.id) {
            baseColumns += [
                'id' : {it.id}, // reserved for tags
            ]
        }
        
        baseColumns += [
            'chr' : {it.chr },
            'pos' : {it.pos },
            'size' : {it.size()},
            'ref': {it.ref },
            'alt': {it.alt },
            'qual': {it.qual },
            'depth': { Variant v ->
                def dp = v.header.samples.grep { v.sampleDosage(it) }.collect { v.getTotalDepth(it) }.max()
                if(dp)
                    return dp
                else
                    return bams*.value*.coverage(v.chr, v.pos).max()
            },
            'vaf' : {THREE_DIGIT_PRECISION.format(it.vaf) },
            'families' : { v ->
                def fcount = v.pedigrees.count {  ped ->
                    def result = ped.samples.any {
                        v.sampleDosage(it)
                    }
                    return result
                }
                return fcount;
            },
            'gqs' : { Variant v -> v.genoTypes.collect { ((it.GQ == null) || it.GQ.equals('.')) ? null : it.GQ.toDouble()} }
        ]

        consColumns = [
            'gene' : {it['SYMBOL']},
            'cons' : {vep -> vep['Consequence'].split('&').min { VepConsequence.fromTerm(it).ordinal() } },
            'maf'  : this.&findMaxMaf
        ]
    }
    
    void processVariant(Variant v,  List<Regions> vcfRegions, Writer w) {
        
        ++stats.total
                    
        for(target in targets) {
            if(!(v in target)) {
                ++stats.excludeByTarget
                return
            }
            
            if(excludeRegions && (v in excludeRegions)) {
                ++stats.excludeByTarget
                return
            }
        }
                    
        List baseInfo = baseColumns.collect { baseColumns[it.key](v) }
        List<Integer> dosages = exportSamples.collect { v.sampleDosage(it) }
//        if(dosages.every { it==0}) {
//            ++stats.excludeNotPresent
//            if(!opts.f) {
//                println "WARNING: variant at $v.chr:$v.pos $v.ref/$v.alt in VCF but not genotyped as present for any sample"
//            }
//            return
//        }
//                        
        if(opts.diff && dosages.clone().unique().size()==1)  {
            ++stats.excludeByDiff
            updateROC(v)
            return
        }
                        
        List<Integer> refCount = v.getAlleleDepths(0)
        List<Integer> altCount = v.getAlleleDepths(1)
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
        consequences.grep { !it.Consequence.split('&').every { excludedVEPConsequences.contains(it) } }.each { vep ->
                        
            if(excludedVEPConsequences.contains(consColumns['cons'](vep))) {
                return
            }
            excludedByCons = false
                        
            if(hasVEP) {
                if(findMaxMaf(vep)>maxMaf)  {
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
            
            if(opts.norep) {
                String base = fasta.basesAt(v.chr, v.pos, v.pos)
                if(base.toLowerCase() == base) {
                    ++stats.excludeByRepeat
                    return
                }
                
                RepeatMotif rep = fasta.repeatAt(v.chr, v.pos)
                if(rep) {
                    if((rep.motif.size() == 1 && rep.repetitions>4) || (rep.motif.size()>1 && rep.repetitions>2)) {
                        ++stats.excludeByRepeat
                        return
                    }
                }
            }
            
//            if(!filters.every { Eval.x(v, it) }) {
            if(!parsedFilters.every { it(v) }) {
                ++stats.excludeByFilter
                return
            }
                        
            if(lastLines>0)
                w.println ","

            List customInfoValues = customInfos.collect { info ->  
                def value = v.parsedInfo[info.infoField] 
                if(!value)
                    return null
                if(info.type == 'list')
                    value = value.tokenize(',').join(', ')
                else
                if(info.type == 'set')
                    value = value.tokenize(',').unique().join(', ')
                if(info.type == 'min') {
                    if(!(value instanceof List))
                        value = value.tokenize(',')
                    value = value.min()
                }
                return value
            }
            List vepInfo = consColumns.collect { name, func -> func(vep) }
            w.print(groovy.json.JsonOutput.toJson(baseInfo+vepInfo+customInfoValues+dosages + [ [refCount,altCount].transpose() ] ))
            printed = true
            
            if(opts.tsv)
                tsvWriter.writeNext((baseInfo+vepInfo+dosages) as String[])
            
            ++lastLines
        }
        
        if(printed) {
            updateROC(v)
            ++stats.totalIncluded
        }
                    
        lastVariant = v;
        if(excludedByCons)
            ++stats.excludeByCons
        
    }
    
    void updateROC(Variant v) {
        if(!opts.roc) 
            return 
        
        int referenceDosage = v.sampleDosage(referenceSample)
        
        for(sample in rocSamples) {
            int rocDosage = v.sampleDosage(sample)
            if(rocDosage>0)
                sampleROCPoints[sample] << new ROCPoint(variant: v, truePositive:rocDosage == referenceDosage)
        }
        
        if(referenceDosage > 0) {
            ++rocTotalVariants
        }
    }

    private printHeadContent(BufferedWriter w) {
        
//        log.info "Sample map: " + Utils.table(sampleMap*.key, [sampleMap*.value])
        w.println """<!DOCTYPE html>\n<html>
                <head>
                <script type='text/javascript'>
            """
        
            if(opts.id)
                w.println "var showId = true;"
                
            if(opts.bams)
                w.println "var bams = " + JsonOutput.toJson(opts.bams)  + ';'

            w.println "var customInfos = " + JsonOutput.toJson(customInfos) + ";"
            w.println "var htmlFilters = " + JsonOutput.toJson(htmlFilters) + ";"
            
           w.println """
                </script>
            """
        w.println css.collect{"<link rel='stylesheet' href='$it'/>"}.join("\n").stripIndent()
        w.println js.collect{"<script type='text/javascript' src='$it'></script>"}.join("\n").stripIndent()

        w.println "<script type='text/javascript'>"
        
        w.println "var variants = ["
    }

    private void injectVCFLibrary(BufferedWriter w) {
        if(opts.vcfjs) {
            w.println "<script type='text/javascript' src='$opts.vcfjs'></script>"
        }
        else {
            // Embed the main vcf.js
            def vcfjs
            def fileVcfJsPath = "src/main/resources/vcf.js"
            if(new File(fileVcfJsPath).exists())
                vcfjs = new File(fileVcfJsPath).text
            else
                vcfjs = this.class.classLoader.getResourceAsStream("vcf.js").text

            w.println "<script type='text/javascript'>\n$vcfjs\n</script>"
        }
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

    @CompileStatic
    List<VCF> loadVCFs() {
        
        List<String> chrs = opts['chrs'] ? (List<String>)opts['chrs'] : null
        
        List<String> preFilters = []
        if(opts['prefilters']) {
            preFilters = (List<String>)opts['prefilters']
            
            log.info "Will apply pre-filters: " + preFilters.join('; ')

            parsedPreFilters = (List<Closure>)preFilters.collect {
                new GroovyShell().evaluate( "{ x ->\n\n$it\n}")
            }
        }
        
        List<String> vcfPaths = (List<String>)opts['is']
        
        List<VCF> vcfs = []
        
        for(String vcfPath in vcfPaths) {
            
            log.info "Read $vcfPath ..."
            
            VCF vcf = new VCF(vcfPath)
            List applicableExportSamples = exportSamples.grep { 
                it in vcf.samples
            }
            
            log.info "The following export samples are applicable to $vcfPath: ${applicableExportSamples.join(',')}"
            
            vcf = VCF.parse(new File(vcfPath), null, samples:applicableExportSamples?:null) { Variant  v ->
                
                if(chrs != null && !(v.chr in chrs))
                    return false
                    
               if(this.targets.any { Regions r -> !r.overlaps(v) }) {
                    ++stats.excludeByTarget
                   return false
               }

                if(!parsedPreFilters.every { it(v) }) {
                    ++stats.excludeByPreFilter
                    return false
                }
            }
            
            log.info "Retained ${vcf.size()} variants from $vcfPath"
            vcfs.add(vcf)
        }
        return vcfs
    }
    
    private void printStats(VCFSummaryStats stats) {
        println " Summary ".center(80,"=")
        
        Writer statsWriter = null
        if(opts.stats)
            statsWriter = new File(opts.stats).newWriter()
            
        try {    
            ["total", "excludeByPreFilter","excludeByDiff", "excludeByCons", "excludeByMaf", "excludeComplex", "excludeByMasked", "excludeByFilter", "excludeByTarget", "excludeByRepeat", "excludeNotPresent","totalIncluded"].each { prop ->
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
        List accumulatedSamples = vcfs[0].samples?.clone()?:['NA']
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
    
    
    
