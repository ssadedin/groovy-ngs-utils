package gngs

import groovy.transform.CompileStatic

/**
 * Models VEP Consequence type as sourced from <a href="https://asia.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences">https://asia.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences</a>
 * (Ensembl release 103 - February 2021).
 *
 * Also merges consequences initiator_codon_variant, non_coding_exon_variant and
 * nc_transcript_variant that can only be found in deprecated {@link VEPConsequences} for backwards
 * compatibility.
 */
@CompileStatic
enum VepConsequence {

    // @formatter:off
    TRANSCRIPT_ABLATION('transcript_ablation'                              , VepImpact.HIGH,     ConsequenceType.UNKNOWN   , 'SO:0001893', '#ff0000', 'Transcript ablation'               , 'A feature ablation whereby the deleted region includes a transcript feature')                                                                                      ,
    SPLICE_ACCEPTOR_VARIANT('splice_acceptor_variant'                      , VepImpact.HIGH,     ConsequenceType.SPLICE    , 'SO:0001574', '#ff581a', 'Splice acceptor variant'           , 'A splice variant that changes the 2 base region at the 3\' end of an intron')                                                                                      ,
    SPLICE_DONOR_VARIANT('splice_donor_variant'                            , VepImpact.HIGH,     ConsequenceType.SPLICE    , 'SO:0001575', '#ff581a', 'Splice donor variant'              , 'A splice variant that changes the 2 base region at the 5\' end of an intron')                                                                                      ,
    STOP_GAINED('stop_gained'                                              , VepImpact.HIGH,     ConsequenceType.NONSENSE  , 'SO:0001587', '#ff0000', 'Stop gained'                       , 'A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript')                       ,
    FRAMESHIFT_VARIANT('frameshift_variant'                                , VepImpact.HIGH,     ConsequenceType.FRAMESHIFT, 'SO:0001589', '#9400d3', 'Frameshift variant'                , 'A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three'),
    STOP_LOST('stop_lost'                                                  , VepImpact.HIGH,     ConsequenceType.NONSENSE  , 'SO:0001578', '#ff0000', 'Stop lost'                         , 'A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript')                                       ,
    START_LOST('start_lost'                                                , VepImpact.HIGH,     ConsequenceType.UNKNOWN   , 'SO:0002012', '#ffd700', 'Start lost'                        , 'A codon variant that changes at least one base of the canonical start codon')                                                                                      ,
    INITIATOR_CODON_VARIANT('initiator_codon_variant'                      , VepImpact.HIGH,     ConsequenceType.SILENT    , 'SO:0001582', '#ffd700', 'Initiator codon variant'           , 'A codon variant that changes at least one base of the first codon of a transcript')                                                                                ,
    TRANSCRIPT_AMPLIFICATION('transcript_amplification'                    , VepImpact.HIGH,     ConsequenceType.UNKNOWN   , 'SO:0001889', '#ff69b4', 'Transcript amplification'          , 'A feature amplification of a region containing a transcript')                                                                                                      ,
    INFRAME_INSERTION('inframe_insertion'                                  , VepImpact.MODERATE, ConsequenceType.MISSENSE  , 'SO:0001821', '#ff69b4', 'Inframe insertion'                 , 'An inframe non synonymous variant that inserts bases into in the coding sequence')                                                                                 ,
    INFRAME_DELETION('inframe_deletion'                                    , VepImpact.MODERATE, ConsequenceType.MISSENSE  , 'SO:0001822', '#ff69b4', 'Inframe deletion'                  , 'An inframe non synonymous variant that deletes bases from the coding sequence')                                                                                    ,
    MISSENSE_VARIANT('missense_variant'                                    , VepImpact.MODERATE, ConsequenceType.MISSENSE  , 'SO:0001583', '#ffd700', 'Missense variant'                  , 'A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved')                               ,
    PROTEIN_ALTERING_VARIANT('protein_altering_variant'                    , VepImpact.MODERATE, ConsequenceType.UNKNOWN   , 'SO:0001818', '#ff0080', 'Protein altering variant'          , 'A sequence_variant which is predicted to change the protein encoded in the coding sequence')                                                                       ,
    SPLICE_REGION_VARIANT('splice_region_variant'                          , VepImpact.LOW,      ConsequenceType.SPLICE    , 'SO:0001630', '#ff7f50', 'Splice region variant'             , 'A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron')           ,
    INCOMPLETE_TERMINAL_CODON_VARIANT('incomplete_terminal_codon_variant'  , VepImpact.LOW,      ConsequenceType.SILENT    , 'SO:0001626', '#ff00ff', 'Incomplete terminal codon variant' , 'A sequence variant where at least one base of the final codon of an incompletely annotated transcript is changed')                                                 ,
    START_RETAINED_VARIANT('start_retained_variant'                        , VepImpact.LOW,      ConsequenceType.UNKNOWN   , 'SO:0002019', '#76ee00', 'Start retained variant'            , 'A sequence variant where at least one base in the start codon is changed, but the start remains')                                                                  ,
    STOP_RETAINED_VARIANT('stop_retained_variant'                          , VepImpact.LOW,      ConsequenceType.SILENT    , 'SO:0001567', '#76ee00', 'Stop retained variant'             , 'A sequence variant where at least one base in the terminator codon is changed, but the terminator remains')                                                        ,
    SYNONYMOUS_VARIANT('synonymous_variant'                                , VepImpact.LOW,      ConsequenceType.SILENT    , 'SO:0001819', '#76ee00', 'Synonymous variant'                , 'A sequence variant where there is no resulting change to the encoded amino acid')                                                                                  ,
    CODING_SEQUENCE_VARIANT('coding_sequence_variant'                      , VepImpact.MODIFIER, ConsequenceType.MISSENSE  , 'SO:0001580', '#458b00', 'Coding sequence variant'           , 'A sequence variant that changes the coding sequence')                                                                                                              ,
    MATURE_MIRNA_VARIANT('mature_miRNA_variant'                            , VepImpact.MODIFIER, ConsequenceType.SILENT    , 'SO:0001620', '#458b00', 'Mature miRNA variant'              , 'A transcript variant located with the sequence of the mature miRNA')                                                                                               ,
    FIVE_PRIME_UTR_VARIANT('5_prime_UTR_variant'                           , VepImpact.MODIFIER, ConsequenceType.SILENT    , 'SO:0001623', '#7ac5cd', '5 prime UTR variant'               , 'A UTR variant of the 5\' UTR')                                                                                                                                     ,
    THREE_PRIME_UTR_VARIANT('3_prime_UTR_variant'                          , VepImpact.MODIFIER, ConsequenceType.SILENT    , 'SO:0001624', '#7ac5cd', '3 prime UTR variant'               , 'A UTR variant of the 3\' UTR')                                                                                                                                     ,
    NON_CODING_TRANSCRIPT_EXON_VARIANT('non_coding_transcript_exon_variant', VepImpact.MODIFIER, ConsequenceType.SILENT    , 'SO:0001792', '#32cd32', 'Non coding transcript exon variant', 'A sequence variant that changes non-coding exon sequence in a non-coding transcript')                                                                              ,
    NON_CODING_EXON_VARIANT('non_coding_exon_variant'                      , VepImpact.MODIFIER, ConsequenceType.SILENT    , 'SO:0001792', '#32cd32', 'Non coding exon variant'           , 'A sequence variant that changes non-coding exon sequence in a non-coding transcript')                                                                              ,
    INTRON_VARIANT('intron_variant'                                        , VepImpact.MODIFIER, ConsequenceType.SILENT    , 'SO:0001627', '#02599c', 'Intron variant'                    , 'A transcript variant occurring within an intron')                                                                                                                  ,
    NMD_TRANSCRIPT_VARIANT('NMD_transcript_variant'                        , VepImpact.MODIFIER, ConsequenceType.UNKNOWN   , 'SO:0001621', '#ff4500', 'NMD transcript variant'            , 'A variant in a transcript that is the target of NMD')                                                                                                              ,
    NON_CODING_TRANSCRIPT_VARIANT('non_coding_transcript_variant'          , VepImpact.MODIFIER, ConsequenceType.UNKNOWN   , 'SO:0001619', '#32cd32', 'Non coding transcript variant'     , 'A transcript variant of a non coding RNA gene')                                                                                                                    ,
    NC_TRANSCRIPT_VARIANT('nc_transcript_variant'                          , VepImpact.MODIFIER, ConsequenceType.SILENT    , 'SO:0001619', '#32cd32', 'Non coding transcript variant'     , 'A transcript variant of a non coding RNA gene')                                                                                                                    ,
    UPSTREAM_GENE_VARIANT('upstream_gene_variant'                          , VepImpact.MODIFIER, ConsequenceType.SILENT    , 'SO:0001631', '#a2b5cd', 'Upstream gene variant'             , 'A sequence variant located 5\' of a gene')                                                                                                                         ,
    DOWNSTREAM_GENE_VARIANT('downstream_gene_variant'                      , VepImpact.MODIFIER, ConsequenceType.SILENT    , 'SO:0001632', '#a2b5cd', 'Downstream gene variant'           , 'A sequence variant located 3\' of a gene')                                                                                                                         ,
    TFBS_ABLATION('TFBS_ablation'                                          , VepImpact.MODIFIER, ConsequenceType.UNKNOWN   , 'SO:0001895', '#a52a2a', 'TFBS ablation'                     , 'A feature ablation whereby the deleted region includes a transcription factor binding site')                                                                       ,
    TFBS_AMPLIFICATION('TFBS_amplification'                                , VepImpact.MODIFIER, ConsequenceType.UNKNOWN   , 'SO:0001892', '#a52a2a', 'TFBS amplification'                , 'A feature amplification of a region containing a transcription factor binding site')                                                                               ,
    TF_BINDING_SITE_VARIANT('TF_binding_site_variant'                      , VepImpact.MODIFIER, ConsequenceType.UNKNOWN   , 'SO:0001782', '#a52a2a', 'TF binding site variant'           , 'A sequence variant located within a transcription factor binding site')                                                                                            ,
    REGULATORY_REGION_ABLATION('regulatory_region_ablation'                , VepImpact.MODERATE, ConsequenceType.UNKNOWN   , 'SO:0001894', '#a52a2a', 'Regulatory region ablation'        , 'A feature ablation whereby the deleted region includes a regulatory region')                                                                                       ,
    REGULATORY_REGION_AMPLIFICATION('regulatory_region_amplification'      , VepImpact.MODIFIER, ConsequenceType.UNKNOWN   , 'SO:0001891', '#a52a2a', 'Regulatory region amplification'   , 'A feature amplification of a region containing a regulatory region')                                                                                               ,
    FEATURE_ELONGATION('feature_elongation'                                , VepImpact.MODIFIER, ConsequenceType.UNKNOWN   , 'SO:0001907', '#7f7f7f', 'Feature elongation'                , 'A sequence variant that causes the extension of a genomic feature, with regard to the reference sequence')                                                         ,
    REGULATORY_REGION_VARIANT('regulatory_region_variant'                  , VepImpact.MODIFIER, ConsequenceType.UNKNOWN   , 'SO:0001566', '#a52a2a', 'Regulatory region variant'         , 'A sequence variant located within a regulatory region')                                                                                                            ,
    FEATURE_TRUNCATION('feature_truncation'                                , VepImpact.MODIFIER, ConsequenceType.UNKNOWN   , 'SO:0001906', '#7f7f7f', 'Feature truncation'                , 'A sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence')                                                         ,
    INTERGENIC_VARIANT('intergenic_variant'                                , VepImpact.MODIFIER, ConsequenceType.SILENT    , 'SO:0001628', '#636363', 'Intergenic variant'                , 'A sequence variant located in the intergenic region, between genes')                                                                                               ,
    UNKNOWN('unknown'                                                      , VepImpact.UNKNOWN , ConsequenceType.UNKNOWN   , 'SO:UNKNOWN', '#636363', 'Unknown'                           , 'Unknown')
    ;
    // @formatter:on

    String term
    VepImpact impact
    ConsequenceType type
    String soAccession  // Sequence Ontology accession number
    String ensemblColourHex
    String displayTerm
    String description

    VepConsequence(String term, VepImpact impact, ConsequenceType type, String soAccession, String ensemblColourHex, String displayTerm, String description) {
        this.term = term
        this.impact = impact
        this.type = type
        this.soAccession = soAccession
        this.ensemblColourHex = ensemblColourHex
        this.displayTerm = displayTerm
        this.description = description
    }

    String getShortName() {
        term.replaceAll('_variant', '')
    }

    int getRank() {
        ordinal() + 1
    }

    int getSeverity() {
        values().size() - ordinal()
    }

    Map<String, Object> toMap() {
        Collections.unmodifiableMap([
                id              : name(),
                rank            : rank,
                severity        : severity,
                name            : name(),
                shortName       : shortName,
                term            : term,
                impact          : impact.toMap(),
                type            : type,
                soAccession     : soAccession,
                ensemblColourHex: ensemblColourHex,
                displayTerm     : displayTerm,
                description     : description,
        ])
    }

    static int severityOfConsequence(VepConsequence consequence) {
        int index = consequence.ordinal()
        if (index < 0)
            index = values().size()
        values().size() - index
    }

    static int severityOf(String consequence) {
        if (consequence == null) {
            return severityOfConsequence(UNKNOWN)
        } else {
            return severityOfConsequence(fromTerm(consequence))
        }
    }

    static knownValues() {
        values() - [UNKNOWN]
    }

    static VepConsequence fromTerm(String term) {
        values().find { it.term == term } ?: UNKNOWN
    }
}

enum ConsequenceType {

    SILENT,
    MISSENSE,
    SPLICE,
    FRAMESHIFT,
    NONSENSE,
    UNKNOWN
    ;

}

enum VepImpact {
    HIGH('The variant is assumed to have high (disruptive) impact in the protein, probably causing protein truncation, loss of function or triggering nonsense mediated decay.'),
    MODERATE('A non-disruptive variant that might change protein effectiveness.'),
    LOW('Assumed to be mostly harmless or unlikely to change protein behaviour.'),
    MODIFIER('Usually non-coding variants or variants affecting non-coding genes, where predictions are difficult or there is no evidence of impact.'),
    UNKNOWN('Unknown')
    ;

    VepImpact(String description) {
        this.description = description
    }

    String description

    Map<String, Object> toMap() {
        return Collections.unmodifiableMap([
                id         : name(),
                name       : name(),
                description: description,
        ])
    }
}
