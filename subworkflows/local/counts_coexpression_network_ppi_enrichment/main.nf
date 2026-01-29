//
// Gene co-expression network analysis, PPI network construction, and functional enrichment
// Performs: WGCNA co-expression network analysis, STRING PPI network analysis, QC visualization, hub gene heatmaps, GO enrichment, KEGG enrichment
//

include { WGCNA_ANALYSIS } from '../../../modules/local/wgcna_analysis/main.nf'
include { QC_VISUALIZATION } from '../../../modules/local/qc_visualization/main.nf'
include { DE_VISUALIZATION } from '../../../modules/local/de_visualization/main.nf'
include { STRING_PPI } from '../../../modules/local/string_ppi/main.nf'
include { GO_ENRICHMENT } from '../../../modules/local/go_enrichment/main.nf'
include { KEGG_ENRICHMENT } from '../../../modules/local/kegg_enrichment/main.nf'

workflow COUNTS_COEXPRESSION_NETWORK_PPI_ENRICHMENT {
    take:
    ch_counts // channel: [ counts_matrix ]
    ch_samplesheet // channel: [ samplesheet ]
    val_min_gs // integer: minimum good samples threshold
    val_min_mm // integer: minimum module membership threshold
    val_top_n_hubgenes // integer: number of top hub genes per module
    val_top_n_processes // integer: number of top processes to visualize
    val_cor_threshold // float: correlation threshold for adjacency
    val_pval_threshold // float: p-value threshold for significance
    val_tomtype // string: TOM type (none/signed/signed 2/unsigned/unsigned 2/signed Nowick/signed Nowick 2)
    val_networktype // string: network type (signed/unsigned/signed hybrid)
    val_deepsplit // integer: deep split level for dendrogram cutting
    val_reassignthreshold // float: reassignment threshold for module merging
    val_mergecutheight // float: merge cut height for dendrogram
    val_minmodulesize // integer: minimum module size
    val_species_name // string: species name for PPI/enrichment
    val_score_threshold // integer: STRING PPI score threshold
    val_enrichment_method // string: enrichment method(s) - "GO", "KEGG", or "GO,KEGG"

    main:

    ch_versions = channel.empty()

    // 1. WGCNA co-expression network analysis
    WGCNA_ANALYSIS(
        ch_counts,
        ch_samplesheet,
        val_min_gs,
        val_min_mm,
        val_top_n_hubgenes,
        val_cor_threshold,
        val_pval_threshold,
        val_tomtype,
        val_networktype,
        val_deepsplit,
        val_reassignthreshold,
        val_mergecutheight,
        val_minmodulesize,
    )
    ch_versions = ch_versions.mix(WGCNA_ANALYSIS.out.versions)

    // 2. STRING protein-protein interaction network analysis
    STRING_PPI(
        WGCNA_ANALYSIS.out.modules_hubgenes,
        val_species_name,
        val_score_threshold,
    )
    ch_versions = ch_versions.mix(STRING_PPI.out.versions)

    // 3. Quality control visualization
    QC_VISUALIZATION(
        [],
        ch_samplesheet,
        WGCNA_ANALYSIS.out.unfiltered_log2counts,
        WGCNA_ANALYSIS.out.filtered_log2counts,
        WGCNA_ANALYSIS.out.normalized_log2counts,
        WGCNA_ANALYSIS.out.vsd_object,
    )
    ch_versions = ch_versions.mix(QC_VISUALIZATION.out.versions)

    // 4. Hub gene heatmap visualization
    DE_VISUALIZATION(
        [],
        ch_samplesheet,
        [],
        0.05,
        1,
        WGCNA_ANALYSIS.out.modules_hubgenes,
    )
    ch_versions = ch_versions.mix(DE_VISUALIZATION.out.versions)

    ch_go_enrichment = channel.empty()
    ch_kegg_enrichment = channel.empty()

    // 5. GO enrichment analysis (optional)
    if ("GO" in val_enrichment_method) {
        GO_ENRICHMENT(
            [],
            val_pval_threshold,
            [],
            val_species_name,
            WGCNA_ANALYSIS.out.modules_hubgenes,
            [],
            [],
            [],
            val_top_n_processes,
        )
        ch_go_enrichment = GO_ENRICHMENT.out.go_results
        ch_versions = ch_versions.mix(GO_ENRICHMENT.out.versions)
    }

    // 6. KEGG pathway enrichment analysis (optional)
    if ("KEGG" in val_enrichment_method) {
        KEGG_ENRICHMENT(
            [],
            val_pval_threshold,
            [],
            val_species_name,
            WGCNA_ANALYSIS.out.modules_hubgenes,
            [],
            [],
            [],
            val_top_n_processes,
        )
        ch_kegg_enrichment = KEGG_ENRICHMENT.out.kegg_results
        ch_versions = ch_versions.mix(KEGG_ENRICHMENT.out.versions)
    }

    emit:
    wgcna_results = WGCNA_ANALYSIS.out.wgcna_results // channel: [ wgcna_results ]
    ppi_results = STRING_PPI.out.string_results // channel: [ ppi_results ]
    qc_plots = QC_VISUALIZATION.out.qc_plots // channel: [ qc_plots ]
    heatmap_plots = DE_VISUALIZATION.out.heatmaps // channel: [ heatmaps ]
    go_enrichment = ch_go_enrichment // channel: [ go_results ]
    kegg_enrichment = ch_kegg_enrichment // channel: [ kegg_results ]
    versions = ch_versions // channel: [ versions.yml ]
}
