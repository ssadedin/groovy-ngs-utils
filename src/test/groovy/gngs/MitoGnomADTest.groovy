package gngs

import groovy.json.JsonOutput
import org.junit.Test

class MitoGnomADTest {

    @Test
    void testParseMitoInfo() {
         def info = """
            qual=-10.0000;variant_collapsed=C5790A;hap_defining_variant;pon_mt_trna_prediction=Neutral;pon_ml_probability_of_pathogenicity=0.298306;mitotip_score=2.11748;mitotip_trna_prediction=Likely Benign;base_qual_hist=0,0,0,0,0,0,0,0,0,0;position_hist=0,0,0,0,0,0,0,0,0,0;strand_bias_hist=0,0,0,0,0,0,0,0,0,0;weak_evidence_hist=0,0,0,0,0,0,0,0,0,0;contamination_hist=1,0,0,0,0,0,0,0,0,0;heteroplasmy_below_10_percent_hist=1,0,0,0,0,0,0,0,0,0;excluded_AC=1;AN=56433;AC_hom=37;AC_het=1;hl_hist=0,0,0,0,0,0,0,1,0,37;hl_hist=0,13,10,4,5,4,8,6,14,35667;dp_mean=3033.10;mq_mean=60.0000;tlod_mean=6297.13;AF_hom=0.000655645;AF_het=1.77201e-05;max_hl=1.00000;hap_AN=2680,1537,868,603,34,282,91,14784,701,934,3144,2732,663,2977,4724,5672,126,1,1297,366,7,393,3080,6037,1234,819,546,12,89;hap_AC_het=0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;hap_AC_hom=0,1,0,0,0,0,0,0,0,0,0,1,2,0,0,0,0,0,0,0,0,0,0,33,0,0,0,0,0;hap_AF_hom=0.00000,0.000650618,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.000366032,0.00301659,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00546629,0.00000,0.00000,0.00000,0.00000,0.00000;hap_AF_het=0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.000318066,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000;hap_hl_hist=[[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,0,0,2],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,33],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0]];hap_faf_hom=0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.000535000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00400003,0.00000,0.00000,0.00000,0.00000,0.00000;hapmax_AF_hom=U;hapmax_AF_het=J;faf_hapmax_hom=0.00400003;age_hist_hom_bin_freq=4,5,1,2,9,4,2,4,1,0;age_hist_hom_n_smaller=1;age_hist_hom_n_larger=0;age_hist_het_bin_freq=0,0,0,0,0,0,0,1,0,0;age_hist_het_n_smaller=0;age_hist_het_n_larger=0;dp_hist_all_n_larger=42621;dp_hist_alt_n_larger=15;dp_hist_all_bin_freq=0,0,41,442,1450,2001,2088,2267,2618,2905;dp_hist_alt_bin_freq=0,0,0,0,4,5,7,2,4,1;vep_csq=A|non_coding_transcript_exon_variant|MODIFIER|MT-TC|ENSG00000210140|Transcript|ENST00000387405|Mt_tRNA|1/1||ENST00000387405.1:n.37G>T||37|||||1||-1|SNV||HGNC|HGNC:7477|YES||||||||||||||||||||
        """.trim()

        def mn = MitoGnomAD.parse(Variant.parseInfoString(info))

        assert mn.ac_het == 1
        assert mn.hl_hist == [0, 13, 10, 4, 5, 4, 8, 6, 14, 35667]

        // hap_AN=2680,1537,868,603,34,282,91,14784,701,934,3144,2732,663,2977,4724,5672,126,1,1297,366,7,393,3080,6037,1234,819,546,12,89
        assert mn.hap_an[-1] == 89
        
        assert mn.hap_defining_variant == true
        
        assert mn.pon_mt_trna_prediction == gngs.PON_TRNA_Prediction.NEUTRAL
        assert mn.mitotip_trna_prediction == MITOTIP_Prediction.LIKELY_BENIGN
        assert mn.pon_ml_probability_of_pathogenicity > 0.2 && mn.pon_ml_probability_of_pathogenicity < 0.4
        
        assert mn.common_low_heteroplasmy == false
        
        // age_hist_hom_bin_freq=4,5,1,2,9,4,2,4,1,0;age_hist_hom_n_smaller=1;age_hist_hom_n_larger=0
        assert mn.age_hist_hom.smaller == 1
        assert mn.age_hist_hom.larger == 0
        assert mn.age_hist_hom.bins[1] == 5
        
        println JsonOutput.prettyPrint(JsonOutput.toJson(mn))
        
    }
}
