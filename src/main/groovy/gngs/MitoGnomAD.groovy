package gngs

import groovy.transform.CompileStatic

@CompileStatic
class GnomADHist {
    
    int smaller
    int larger
    
    List<Integer> bins;
    
    static GnomADHist parse(def value, def smaller, def larger) {
        new GnomADHist(bins: MitoGnomAD.intList(value), smaller: smaller as Integer, larger: larger as Integer)
    }
}

enum PON_TRNA_Prediction {
    PATHOGENIC, LIKELY_PATHOGENIC, LIKELY_NEUTRAL, NEUTRAL, VUS
}

enum MITOTIP_Prediction {
    LIKELY_PATHOGENIC, POSSIBLY_PATHOGENIC, POSSIBLY_BENIGN, LIKELY_BENIGN
}

/**
 * Models gnomAD annotations for mitochondrial variants
 * 
 * @author Simon Sadedin
 */
class MitoGnomAD {
    
    public static final List<String> HAPLOGROUPS= ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'HV', 'I', 'J', 'K', 'L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'M', 'N', 'P', 'R', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
    
   
   @CompileStatic
    static List<Integer> intList(def value) {
        String.valueOf(value).tokenize(',')*.toInteger()
    }
    
    @CompileStatic
    static List<Double> doubleList(def value) {
        String.valueOf(value).tokenize(',')*.toDouble()
    }
    
    @CompileStatic
    static List<List<Integer>> listOfHistograms(def value) {
        ((String)value)[1..-2].tokenize('],[').collate(10).collect { it*.toInteger() }
    }
  
    
    String variant_collapsed

    boolean hap_defining_variant
    boolean common_low_heteroplasmy

    GnomADHist base_qual_hist;
    GnomADHist position_hist;
    GnomADHist strand_bias_hist;
    GnomADHist weak_evidence_hist;
    GnomADHist contamination_hist;
    GnomADHist heteroplasmy_below_10_percent_hist;
    int excluded_AC
    int an
    int ac_hom
    int ac_het
    GnomADHist hl_hist;
    double dp_mean;
    double mq_mean;
    double tlod_mean;
    double af_hom;
    double af_het;

    /**@
     * Maximum heteroplasmy level observed among all samples for that variant
     */
    double max_hl;

    List<Integer> hap_an

    List<Integer> hap_ac_het;
    List<Integer> hap_ac_hom;
    List<Double> hap_af_hom;
    List<Double> hap_af_het;
    List<List<Integer>> hap_hl_hist;
    List<Double> hap_faf_hom;
    String hapmax_af_hom
    double faf_hapmax_hom;
    GnomADHist age_hist_hom;
    GnomADHist age_hist_het;
    GnomADHist dp_hist_all;
    GnomADHist dp_hist_alt;
    
    // Optional fields: PON TRNA predictions
    PON_TRNA_Prediction pon_mt_trna_prediction
    Double pon_ml_probability_of_pathogenicity
    
    // Optional fields: MITOTIP scores
    Double mitotip_score
    MITOTIP_Prediction mitotip_trna_prediction 
   
    @CompileStatic
    static MitoGnomAD parse(Map<String,Object> info) {
        MitoGnomAD result = new MitoGnomAD(
            an : info.AN as Integer,
            ac_hom : info.AC_hom as Integer,
            ac_het : info.AC_het as Integer,
            tlod_mean : info.tlod_mean as Double,
            af_hom : info.AF_hom as Double,
            af_het : info.AF_het as Double,
            max_hl : info.max_hl as Double,
            hap_an : intList(info.hap_AN),
            hap_ac_het : intList(info.hap_AC_het),
            hap_ac_hom : intList(info.hap_AC_hom),
            hap_af_hom : doubleList(info.hap_AF_hom),
            hap_af_het : doubleList(info.hap_AF_het),
            hap_hl_hist : listOfHistograms(info.hap_hl_hist),
            hap_faf_hom : doubleList(info.hap_faf_hom),
            hapmax_af_hom : (String)info.hapmax_AF_hom,
            faf_hapmax_hom : info.faf_hapmax_hom as Double,
            age_hist_hom : GnomADHist.parse(info.age_hist_hom_bin_freq, info.age_hist_hom_n_smaller, info.age_hist_hom_n_larger),
            age_hist_het : GnomADHist.parse(info.age_hist_het_bin_freq, info.age_hist_het_n_smaller, info.age_hist_het_n_larger),
            dp_hist_all : GnomADHist.parse(info.dp_hist_all_bin_freq, 0, info.dp_hist_all_n_larger),
            dp_hist_alt : GnomADHist.parse(info.dp_hist_alt_bin_freq, 0, info.dp_hist_alt_n_larger),
            hap_defining_variant : info.containsKey('hap_defining_variant'),
            pon_mt_trna_prediction : info.pon_mt_trna_prediction ? 
                PON_TRNA_Prediction.valueOf(((String)info.pon_mt_trna_prediction).toUpperCase().replaceAll(' ','_'))
                :
                null,
            pon_ml_probability_of_pathogenicity : info.pon_ml_probability_of_pathogenicity ?
                info.pon_ml_probability_of_pathogenicity as Double : null,
                
            mitotip_trna_prediction : info.mitotip_trna_prediction ?
                MITOTIP_Prediction.valueOf(((String)info.mitotip_trna_prediction).toUpperCase().replaceAll(' ','_'))
                : null,
            common_low_heteroplasmy: info.containsKey('common_low_heteroplasmy')
        )
    }
}
