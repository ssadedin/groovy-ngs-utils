import groovy.xml.StreamingMarkupBuilder

int count = 1

Cli cli = new Cli()
cli.with {
    i 'vcf File', args:1, required:true
    p 'PED file describing relationships between the samples in the data', args:1
    o 'Output HTML file', args:1, required:true
    f 'Comma separated list of families to export', args:1
    maxMaf 'Filter out variants above this MAF', args:1
}

opts = cli.parse(args)

Pedigrees pedigrees = null

if(opts.p) {
    pedigrees = Pedigrees.parse(opts.p)
}
else {
    pedigrees = Pedigrees.fromSingletons(new VCF(opts.i).samples) 
}

def exportSamples = pedigrees.families.values().collect { it.individuals*.id }.flatten().grep { it in new VCF(opts.i).samples };
def exportFamilies = pedigrees.families.keySet()
if(opts.f) {
    exportFamilies = opts.f.split(",").collect { it.trim() }
    exportSamples = exportSamples.grep { s -> exportFamilies.any { f-> pedigrees.families[f].samples.contains(s)  }}
    pedigrees.families.keySet().grep { !exportFamilies.contains(it) }.each { pedigrees.removeFamily(it) }
}
println "Samples to export: $exportSamples" 

float MAF_THRESHOLD=0.05
if(opts.maxMaf) 
    MAF_THRESHOLD=opts.maxMaf.toFloat()

def js = [
    "http://ajax.aspnetcdn.com/ajax/jQuery/jquery-2.1.0.min.js",
    "http://cdn.datatables.net/1.10.0/js/jquery.dataTables.min.js",
    "http://ajax.googleapis.com/ajax/libs/jqueryui/1.10.4/jquery-ui.min.js",
    "http://cdnjs.cloudflare.com/ajax/libs/jquery-layout/1.3.0-rc-30.79/jquery.layout.min.js",
    "vcf.js"
]

def css = [
    "http://cdn.datatables.net/1.10.0/css/jquery.dataTables.css",
    "//ajax.googleapis.com/ajax/libs/jqueryui/1.10.4/themes/smoothness/jquery-ui.css",
    "http://cdnjs.cloudflare.com/ajax/libs/jquery-layout/1.3.0-rc-30.79/layout-default.css"
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

def json(obj) {
    groovy.json.JsonOutput.toJson(obj)
}

findMaxMaf = { vep -> 
    [vep.EA_MAF, vep.ASN_MAF, vep.EUR_MAF].collect{it?it.split('&'):[]}.flatten().collect { it.toFloat()}.max() ?: 0.0f
}

def baseColumns = new LinkedHashMap()
baseColumns += [
    'chr' : {it.chr },
    'pos' : {it.pos },
    'ref': {it.ref },
    'alt': {it.alt },
    'qual': {it.qual },
    'depth': {it.info.DP },
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

def vepColumns = [ 
    'gene' : {it['SYMBOL']},
    'cons' : {vep -> vep['Consequence'].split('&').min { VEP_PRIORITY.indexOf(it) } },
    'maf'  : findMaxMaf
]

new File(opts.o).withWriter { w ->
    w.println """
    <html>
        <head>
    """
    w.println css.collect{"<link rel='stylesheet' href='$it'/>"}.join("\n").stripIndent()
    w.println js.collect{"<script type='text/javascript' src='$it'></script>"}.join("\n").stripIndent()
    
    w.println "<script type='text/javascript'>"
   
    w.println "var variants = ["
    int i=0;
    int lastLines = 0
    Variant last = null;
    VCF.parse(opts.i, pedigrees) { v->
        List baseInfo = baseColumns.collect { baseColumns[it.key](v) }
        List dosages = exportSamples.collect { v.sampleDosage(it) }
        if(dosages.every { it==0})
            return
        
        if(i++>0 && lastLines > 0)
            w.println ","
        
        List<Map> veps = v.vepInfo
        lastLines = 0
        veps.grep { !it.Consequence.split('&').every { EXCLUDE_VEP.contains(it) } }.each { vep ->
            
            if(EXCLUDE_VEP.contains(vepColumns['cons'](vep)))
                return
                
            if(findMaxMaf(vep)>MAF_THRESHOLD) 
                return
            
            if(lastLines>0)
                w.println ","
                
            List vepInfo = vepColumns.collect { name, func -> func(vep) }
            w.print(groovy.json.JsonOutput.toJson(baseInfo+vepInfo+dosages))
            ++lastLines
        }
        last = v;
    }
    w.println """];"""
    
    w.println "var columnNames = ${json(baseColumns*.key + vepColumns*.key + exportSamples)};"
    
    w.println """
    var samples = ${json(exportSamples)};

    var pedigrees = ${pedigrees.toJson()};

    var variantTable = null;
    \$(document).ready(function() {
        variantTable = \$.VariantTable('variantTable', samples, variants);
    });
    </script>
    <style type='text/css'>
    td.vcfcol { text-align: center; }
    </style>
    """;
    
   
    w.println """
    </head>
    <body>
    <p>Loading ...</p>
        <div class="ui-layout-north">
            <h1>VCF File ${new File(opts.i).name}</h1>
            <span id=filterOuter><span id=filters> </span></span>
            <div id=filterHelp style='display:none'> </div>
        </div>
        <div class="ui-layout-center">
            <div id=tableHolder>
            </div>
        </div>
        <div class="ui-layout-south">Report created ${new Date()}</div>
        <div class="ui-layout-east" id="ui-layout-east"></div>
    </body>
    </html>
    """
}