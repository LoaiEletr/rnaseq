//
// Isoform switch analysis with functional enrichment
// Performs: IsoformSwitchAnalyzeR analysis for differential isoform usage/alternative splicing,
// identification of switching isoforms, GO enrichment, KEGG enrichment
//

include { ISOFORMSWITCHANALYZER } from '../../../modules/local/isoformswitchanalyzer/main.nf'
include { GO_ENRICHMENT } from '../../../modules/local/go_enrichment/main.nf'
include { KEGG_ENRICHMENT } from '../../../modules/local/kegg_enrichment/main.nf'

workflow COUNTS_ISOFORM_SWITCH_ENRICHMENT {
    take:
    ch_counts // channel: [ quant_results ]
    ch_samplesheet // channel: [ samplesheet ]
    ch_gtf // channel: [ gtf ]
    ch_transcript // channel: [ transcript_fasta ]
    val_quant_type // string: quantification type - "kallisto" or "salmon"
    val_method // string: analysis method - "DIU", "AS", or combinations
    val_pvalue_threshold // float: p-value threshold for significance
    val_dif_cutoff // float: log fold change threshold
    val_ntop_isoforms // integer: number of top isoforms to output
    val_species_name // string: species name for enrichment analysis
    val_ntop_processes // integer: number of top processes to output
    val_enrichment_method // string: enrichment method(s) - "GO", "KEGG", or combinations

    main:
    ch_versions = channel.empty()
    ch_go_enrichment = channel.empty()
    ch_kegg_enrichment = channel.empty()

    // 1. Isoform switch analysis
    ISOFORMSWITCHANALYZER(
        ch_counts,
        ch_samplesheet,
        val_quant_type,
        ch_gtf,
        ch_transcript,
        val_method,
        val_dif_cutoff,
        val_pvalue_threshold,
        val_ntop_isoforms,
    )
    ch_versions = ch_versions.mix(ISOFORMSWITCHANALYZER.out.versions)

    // 2. GO enrichment analysis (optional)
    if ("GO" in val_enrichment_method) {
        GO_ENRICHMENT(
            [],
            val_pvalue_threshold,
            [],
            val_species_name,
            [],
            [],
            ISOFORMSWITCHANALYZER.out.significant_geneids.ifEmpty([]),
            [],
            val_ntop_processes,
        )
        ch_go_enrichment = GO_ENRICHMENT.out.go_results
        ch_versions = ch_versions.mix(GO_ENRICHMENT.out.versions)
    }

    // 3. KEGG pathway enrichment analysis (optional)
    if ("KEGG" in val_enrichment_method) {
        KEGG_ENRICHMENT(
            [],
            val_pvalue_threshold,
            [],
            val_species_name,
            [],
            [],
            ISOFORMSWITCHANALYZER.out.significant_geneids.ifEmpty([]),
            [],
            val_ntop_processes,
        )
        ch_kegg_enrichment = KEGG_ENRICHMENT.out.kegg_results
        ch_versions = ch_versions.mix(KEGG_ENRICHMENT.out.versions)
    }

    emit:
    isoform_results = ISOFORMSWITCHANALYZER.out.isoform_results // channel: [ isoform_results ]
    significant_geneids = ISOFORMSWITCHANALYZER.out.significant_geneids // channel: [ significant_geneids ]
    go_enrichment = ch_go_enrichment // channel: [ go_results ]
    kegg_enrichment = ch_kegg_enrichment // channel: [ kegg_results ]
    versions = ch_versions // channel: [ versions.yml ]
}
