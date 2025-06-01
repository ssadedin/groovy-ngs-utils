package gngs

import graxxia.FileLike
import groovy.transform.CompileStatic
import groovy.util.logging.Log



/**
 * A helper class that provides an interface to easily extract
 * information from downloadable ClinVar resources.
 * 
 * @author simon.sadedin
 */
@Log
class ClinVar {
    
    public final static String CITATIONS_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/var_citations.txt"
    
    public final static String VCF_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh##grchVersion##/clinvar.vcf.gz"
    
    HashMap<String,HashMap<String,List<Integer>>> alleleCitations
    
    HashMap<String, List<Variant>> geneToPathogenicVariants
    
    def vcfPath

    def citationsPath
    
    /**
     * TODO
     * 
     * want to enable pathway of Gene => Variants in Gene => Citations for Variants in Gene
     * 
     * SO: need to 
     * 
     * - [DONE] add resource for ClinVar VCF file
     * - [DONE] parse ClinVar VCF file
     *   [DONE] - index variants by gene symbol
     * - [DONE] write unit test for above functions
     * - refactor variant info into ClinvarVariant class with typed fields / accessors?
     * - add function that finds variants in gene with given pathogencity from gene and then
     *   queries their pubmed article ids
     * - add caching so the load is much faster?
     * - go to 
     */
    
    ClinVar(def vcfPath, def citationsPath) {
        this.vcfPath = vcfPath
        this.citationsPath = citationsPath 
    }
    
    final static Set<String> UNCERTAIN_VARIANT_EVIDENCE = [
        'no_assertion_criteria_provided',
        'no_assertion_provided',
        'no_interpretation_for_the_single_variant',
        'reviewed_by_expert_panel'
    ] as Set
    
    @CompileStatic
    ClinVar load(Map options=[:]) {
        
        final boolean filterToPathogenicOnly = (options.pathogenicOnly == true)
        final boolean filterToHighConfidenceOnly = (options.highConfidenceOnly == true)
        
        ClinVar me = this

        FileLike.reader(citationsPath) { reader ->
            def result = new graxxia.TSV(reader).groupBy {
                it['#AlleleID']
            }
            .collectEntries {
              [ it.key, it.value.groupBy { it['citation_source'] }.collectEntries { [it.key, it.value['citation_id']] }  ]
            }
            
            me.alleleCitations = (HashMap<String,HashMap<String,List<Integer>>>)result
        }
        
        this.geneToPathogenicVariants = new LinkedHashMap(20000)
        
//        File cachedPath = new File(vcfPath.toString() + ".ser")
//        if(vcfPath instanceof String) {
//            if(cachedPath.exists()) {
//                log.info "Loading variant index from cache ...."
//                me.geneToPathogenicVariants = (HashMap<String, List<Variant>>)cachedPath.withObjectInputStream { ois -> ois.readObject() }
//            }
//        }
//
        FileLike.reader(vcfPath) { reader ->
            VCF.parse(reader, false) { v ->
                if(filterToPathogenicOnly) {
                    boolean isPathogenic = ((String)v.info.CLNSIG)?.toLowerCase()?.contains('pathogenic')?:false
                    if(!isPathogenic) {
                        return false
                    }
                }
                
                if(filterToHighConfidenceOnly) {
                    if(v.info.CLNREVSTAT in UNCERTAIN_VARIANT_EVIDENCE) {
                        return false
                    }
                }

                // Passes filtering, add to variant index
                String gene = ((String)v.info.GENEINFO)?.tokenize(':')?.getAt(0)
                if(!gene)
                    return false
                if(me.geneToPathogenicVariants.containsKey(gene))
                    me.geneToPathogenicVariants[gene].add(v)
                else
                    me.geneToPathogenicVariants[gene] = [v]
            }
        }

//        if(vcfPath instanceof String && options.cache) {
//            log.info "Writing cache of variant index to $cachedPath"
//            cachedPath.withObjectOutputStream { oos ->  oos << me.geneToPathogenicVariants }
//        }
        
        log.info "Variants loaded from $vcfPath"
        
        return me
    }

    @CompileStatic
    static ClinVar download(String genomeVersion="hg19") {
        GenomeResource resource = new ResourceDownloader(CITATIONS_URL).download(genomeVersion)
        GenomeResource vcfResource = new ResourceDownloader(VCF_URL).download(genomeVersion)
        return new ClinVar(vcfResource.path, resource.path)
    }
}
