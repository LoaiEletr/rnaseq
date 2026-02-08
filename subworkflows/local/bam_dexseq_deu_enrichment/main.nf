//
// DEXSeq differential exon usage analysis with functional enrichment
// Performs: exon-level counting, differential exon usage analysis,
// MA plots and gene plots, GO enrichment, KEGG enrichment
//

include { DEXSEQ_COUNT } from '../../../modules/local/dexseq/count/main.nf'
include { DEXSEQ_DEU } from '../../../modules/local/dexseq/deu/main.nf'
include { GO_ENRICHMENT } from '../../../modules/local/go_enrichment/main.nf'
include { KEGG_ENRICHMENT } from '../../../modules/local/kegg_enrichment/main.nf'

workflow BAM_DEXSEQ_DEU_ENRICHMENT {
    take:
    ch_bam // channel: [ val(meta), bam ]
    ch_gff // channel: [ gff ]
    ch_samplesheet // channel: [ samplesheet ]
    val_alignment_quality // integer: minimum alignment quality score
    val_pval_threshold // float: p-value threshold for significance
    val_logfc_threshold // float: log fold change threshold
    val_top_n // integer: number of top significant exons to output
    val_min_exonlength // integer: minimum exon length to consider
    val_species_name // string: species name for enrichment analysis
    val_ntop_processes // integer: number of top processes to output
    val_enrichment_method // string: enrichment method(s) - "GO", "KEGG", or "GO,KEGG"

    main:
    ch_versions = channel.empty()
    ch_go_enrichment = channel.empty()
    ch_kegg_enrichment = channel.empty()

    // 1. Count exon-level reads from BAM files
    DEXSEQ_COUNT(
        ch_bam,
        ch_gff,
        val_alignment_quality,
    )
    ch_counts = DEXSEQ_COUNT.out.counts
    ch_versions = ch_versions.mix(DEXSEQ_COUNT.out.versions.first())

    // 2. Differential exon usage analysis
    DEXSEQ_DEU(
        ch_counts.map { meta, counts ->
            counts
        }.collect(),
        ch_samplesheet,
        ch_gff,
        val_pval_threshold,
        val_logfc_threshold,
        val_top_n,
        val_min_exonlength,
    )
    ch_versions = ch_versions.mix(DEXSEQ_DEU.out.versions)

    // 3. GO enrichment analysis (optional)
    if ("GO" in val_enrichment_method) {
        GO_ENRICHMENT(
            [],
            val_pval_threshold,
            [],
            val_species_name,
            [],
            [],
            [],
            DEXSEQ_DEU.out.significant_genes,
            val_ntop_processes,
        )
        ch_go_enrichment = GO_ENRICHMENT.out.go_results
        ch_versions = ch_versions.mix(GO_ENRICHMENT.out.versions)
    }

    // 4. KEGG pathway enrichment analysis (optional)
    if ("KEGG" in val_enrichment_method) {
        KEGG_ENRICHMENT(
            [],
            val_pval_threshold,
            [],
            val_species_name,
            [],
            [],
            [],
            DEXSEQ_DEU.out.significant_genes,
            val_ntop_processes,
        )
        ch_kegg_enrichment = KEGG_ENRICHMENT.out.kegg_results
        ch_versions = ch_versions.mix(KEGG_ENRICHMENT.out.versions)
    }

    emit:
    dexseq_table = DEXSEQ_DEU.out.results_table // channel: [ results_table ]
    dexseq_sigtable = DEXSEQ_DEU.out.significant_results // channel: [ significant_results ]
    ma_plot = DEXSEQ_DEU.out.ma_plot // channel: [ ma_plot ]
    deu_plots = DEXSEQ_DEU.out.gene_plots // channel: [ gene_plots ]
    significant_genes = DEXSEQ_DEU.out.significant_genes // channel: [ significant_genes ]
    report = DEXSEQ_DEU.out.report // channel: [ report ]
    go_enrichment = ch_go_enrichment // channel: [ go_results ]
    kegg_enrichment = ch_kegg_enrichment // channel: [ kegg_results ]
    versions = ch_versions // channel: [ versions.yml ]
}
