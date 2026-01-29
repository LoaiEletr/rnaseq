//
// Differential expression analysis pipeline with quality control, visualization, and functional enrichment
// Performs: edgeR/limma (voom) or DESeq2 differential expression analysis, maSigPro timecourse analysis,
// quality control visualizations, DE visualization (volcano/heatmaps), GO enrichment, KEGG enrichment, MSigDB GSEA
//

include { EDGER_VOOM } from '../../../modules/local/edger_voom/main.nf'
include { LIMMA_DE } from '../../../modules/local/limma_de/main.nf'
include { DESEQ2_DE } from '../../../modules/local/deseq2_de/main.nf'
include { MASIGPRO_TIMECOURSE } from '../../../modules/local/masigpro_timecourse/main.nf'
include { QC_VISUALIZATION } from '../../../modules/local/qc_visualization/main.nf'
include { DE_VISUALIZATION } from '../../../modules/local/de_visualization/main.nf'
include { GO_ENRICHMENT } from '../../../modules/local/go_enrichment/main.nf'
include { KEGG_ENRICHMENT } from '../../../modules/local/kegg_enrichment/main.nf'
include { MSIGDB_GSEA } from '../../../modules/local/msigdb_gsea/main.nf'

workflow COUNTS_DIFFEXPR_QC_ENRICHMENT_VISUALIZATION {
    take:
    ch_counts // channel: [ counts_matrix ]
    ch_samplesheet // channel: [ samplesheet ]
    val_pvalue_threshold // float: p-value threshold for significance
    val_logfc_threshold // float: log fold change threshold
    val_ntop_genes // integer: number of top genes to output
    val_diffexpr_method // string: differential expression method - "limma", "deseq2", or "masigpro"
    val_species_name // string: species name for enrichment analysis
    val_nes_threshold // float: normalized enrichment score threshold for GSEA
    val_msigdb_categories // string: MSigDB categories for GSEA (comma-separated)
    val_padj_gsea // float: adjusted p-value threshold for GSEA
    val_ntop_processes // integer: number of top processes to output
    val_rank_method // string: ranking method for GSEA - "t_stat", "logfc", or "signed_significance"
    val_rsq_threshold // float: R-squared threshold for maSigPro
    val_cluster_method // string: clustering method for maSigPro - "hclust", "Mclust", or "kmeans"
    val_enrichment_method // string: enrichment method(s) - "GO", "KEGG", "GSEA", or combinations

    main:

    ch_versions = channel.empty()

    // Initialize empty channels for conditional outputs
    ch_unfiltered_log2counts = channel.empty()
    ch_filtered_log2counts = channel.empty()
    ch_normalized_log2counts = channel.empty()
    ch_voom_object = channel.empty()
    ch_voom_plot = channel.empty()
    ch_vst_object = channel.empty()
    ch_ebfit_rds = channel.empty()
    ch_diffgenes_expression = channel.empty()
    ch_diffgenes_table = channel.empty()
    ch_unfiltered_deg = channel.empty()
    ch_topgenes_rds = channel.empty()
    ch_topgenes_csv = channel.empty()
    ch_timecourse_clusterplots = channel.empty()
    ch_timecourse_clustermedian = channel.empty()
    ch_timecourse_summarystatistics = channel.empty()
    ch_timecourse_clusterassignments = channel.empty()
    ch_timecourse_clusterlist = channel.empty()
    ch_go_enrichment = channel.empty()
    ch_kegg_enrichment = channel.empty()
    ch_gsea_results = channel.empty()

    // 1. Differential expression analysis (choose one method)
    if (val_diffexpr_method == "limma") {
        // edgeR/limma-voom analysis
        EDGER_VOOM(
            ch_counts,
            ch_samplesheet,
        )
        ch_unfiltered_log2counts = EDGER_VOOM.out.unfiltered_log2cpm
        ch_filtered_log2counts = EDGER_VOOM.out.filtered_log2cpm
        ch_normalized_log2counts = EDGER_VOOM.out.normalized_log2cpm
        ch_voom_object = EDGER_VOOM.out.voom_transformed
        ch_voom_plot = EDGER_VOOM.out.voom_plot
        ch_versions = ch_versions.mix(EDGER_VOOM.out.versions)

        LIMMA_DE(
            ch_voom_object,
            ch_samplesheet,
            val_pvalue_threshold,
            val_logfc_threshold,
            val_ntop_genes,
        )
        ch_ebfit_rds = LIMMA_DE.out.ebfit_rds
        ch_diffgenes_expression = LIMMA_DE.out.de_expression_rds
        ch_diffgenes_table = LIMMA_DE.out.de_genes_rds
        ch_unfiltered_deg = LIMMA_DE.out.unfiltered_deg_rds
        ch_topgenes_rds = LIMMA_DE.out.topgenes_rds
        ch_topgenes_csv = LIMMA_DE.out.topgenes_csv
        ch_versions = ch_versions.mix(LIMMA_DE.out.versions)
    }
    else if (val_diffexpr_method == "deseq2") {
        // DESeq2 analysis
        DESEQ2_DE(
            ch_counts,
            ch_samplesheet,
            val_pvalue_threshold,
            val_logfc_threshold,
            val_ntop_genes,
        )
        ch_unfiltered_log2counts = DESEQ2_DE.out.unfiltered_log2counts
        ch_filtered_log2counts = DESEQ2_DE.out.filtered_log2counts
        ch_normalized_log2counts = DESEQ2_DE.out.normalized_log2counts
        ch_vst_object = DESEQ2_DE.out.vst_object
        ch_diffgenes_table = DESEQ2_DE.out.diffgenes_table
        ch_diffgenes_expression = DESEQ2_DE.out.diffgenes_matrix
        ch_unfiltered_deg = DESEQ2_DE.out.unfiltered_deg
        ch_topgenes_rds = DESEQ2_DE.out.topgenes_rds
        ch_topgenes_csv = DESEQ2_DE.out.topgenes_csv
        ch_versions = ch_versions.mix(DESEQ2_DE.out.versions)
    }
    else if (val_diffexpr_method == "masigpro") {
        // maSigPro timecourse analysis
        MASIGPRO_TIMECOURSE(
            ch_counts,
            ch_samplesheet,
            val_pvalue_threshold,
            val_rsq_threshold,
            val_cluster_method,
        )
        ch_unfiltered_log2counts = MASIGPRO_TIMECOURSE.out.unfiltered_log2cpm
        ch_filtered_log2counts = MASIGPRO_TIMECOURSE.out.filtered_log2cpm
        ch_normalized_log2counts = MASIGPRO_TIMECOURSE.out.normalized_log2cpm
        ch_timecourse_clusterplots = MASIGPRO_TIMECOURSE.out.cluster_plots
        ch_timecourse_clustermedian = MASIGPRO_TIMECOURSE.out.cluster_median_expression
        ch_timecourse_summarystatistics = MASIGPRO_TIMECOURSE.out.summary_statistics
        ch_timecourse_clusterassignments = MASIGPRO_TIMECOURSE.out.cluster_assignments
        ch_timecourse_clusterlist = MASIGPRO_TIMECOURSE.out.clusters_list
        ch_diffgenes_expression = MASIGPRO_TIMECOURSE.out.significant_expression
        ch_versions = ch_versions.mix(MASIGPRO_TIMECOURSE.out.versions)
    }
    else {
        error("Invalid diffexpr_method: ${val_diffexpr_method}. Must be one of: 'limma', 'deseq2', or 'masigpro'")
    }

    // 2. Quality control visualization
    QC_VISUALIZATION(
        ch_voom_object.ifEmpty([]),
        ch_samplesheet,
        ch_unfiltered_log2counts,
        ch_filtered_log2counts,
        ch_normalized_log2counts,
        ch_vst_object.ifEmpty([]),
    )
    ch_versions = ch_versions.mix(QC_VISUALIZATION.out.versions)

    // 3. Differential expression visualization
    DE_VISUALIZATION(
        ch_unfiltered_deg.ifEmpty([]),
        ch_samplesheet,
        ch_timecourse_clusterlist.ifEmpty([]),
        val_pvalue_threshold,
        val_logfc_threshold,
        ch_diffgenes_expression,
    )
    ch_versions = ch_versions.mix(DE_VISUALIZATION.out.versions)

    // 4. GO enrichment analysis (optional)
    if ("GO" in val_enrichment_method) {
        GO_ENRICHMENT(
            ch_unfiltered_deg.ifEmpty([]),
            val_pvalue_threshold,
            val_logfc_threshold,
            val_species_name,
            [],
            ch_timecourse_clusterlist.ifEmpty([]),
            [],
            [],
            val_ntop_processes,
        )
        ch_go_enrichment = GO_ENRICHMENT.out.go_results
        ch_versions = ch_versions.mix(GO_ENRICHMENT.out.versions)
    }

    // 5. KEGG pathway enrichment analysis (optional)
    if ("KEGG" in val_enrichment_method) {
        KEGG_ENRICHMENT(
            ch_unfiltered_deg.ifEmpty([]),
            val_pvalue_threshold,
            val_logfc_threshold,
            val_species_name,
            [],
            ch_timecourse_clusterlist.ifEmpty([]),
            [],
            [],
            val_ntop_processes,
        )
        ch_kegg_enrichment = KEGG_ENRICHMENT.out.kegg_results
        ch_versions = ch_versions.mix(KEGG_ENRICHMENT.out.versions)
    }

    // 6. MSigDB gene set enrichment analysis (optional)
    if ("GSEA" in val_enrichment_method) {
        MSIGDB_GSEA(
            ch_unfiltered_deg,
            val_pvalue_threshold,
            val_species_name,
            val_msigdb_categories,
            val_nes_threshold,
            val_padj_gsea,
            val_ntop_processes,
            val_rank_method,
        )
        ch_gsea_results = MSIGDB_GSEA.out.gsea_results
        ch_versions = ch_versions.mix(MSIGDB_GSEA.out.versions)
    }

    emit:
    unfiltered_log2counts = ch_unfiltered_log2counts // channel: [ unfiltered_log2counts ]
    filtered_log2counts = ch_filtered_log2counts // channel: [ filtered_log2counts ]
    normalized_log2counts = ch_normalized_log2counts // channel: [ normalized_log2counts ]
    voom_object = ch_voom_object // channel: [ voom_object ]
    voom_plot = ch_voom_plot // channel: [ voom_plot ]
    ebfit_rds = ch_ebfit_rds // channel: [ ebfit_rds ]
    diffgenes_expression = ch_diffgenes_expression // channel: [ diffgenes_expression ]
    diffgenes_table = ch_diffgenes_table // channel: [ diffgenes_table ]
    unfiltered_deg = ch_unfiltered_deg // channel: [ unfiltered_deg ]
    topgenes_rds = ch_topgenes_rds // channel: [ topgenes_rds ]
    topgenes_csv = ch_topgenes_csv // channel: [ topgenes_csv ]
    timecourse_clusterplots = ch_timecourse_clusterplots // channel: [ cluster_plots ]
    timecourse_clustermedian = ch_timecourse_clustermedian // channel: [ cluster_median_expression ]
    timecourse_summarystatistics = ch_timecourse_summarystatistics // channel: [ summary_statistics ]
    timecourse_clusterassignments = ch_timecourse_clusterassignments // channel: [ cluster_assignments ]
    timecourse_clusterlist = ch_timecourse_clusterlist // channel: [ clusters_list ]
    qc_plots = QC_VISUALIZATION.out.qc_plots // channel: [ qc_plots ]
    volcano_plot = DE_VISUALIZATION.out.volcano_plot // channel: [ volcano_plot ]
    heatmap_plots = DE_VISUALIZATION.out.heatmaps // channel: [ heatmaps ]
    go_enrichment = ch_go_enrichment // channel: [ go_results ]
    kegg_enrichment = ch_kegg_enrichment // channel: [ kegg_results ]
    gsea_results = ch_gsea_results // channel: [ gsea_results ]
    versions = ch_versions // channel: [ versions.yml ]
}
