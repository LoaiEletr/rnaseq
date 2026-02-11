#!/usr/bin/env nextflow

process WGCNA_ANALYSIS {
    tag "${count_matrix_rds}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/30/30159b6f9d77a89104519ad9f53a798367d227065a206d47afba0eecfbeef78c/data'
        : 'community.wave.seqera.io/library/bioconductor-deseq2_r-wgcna_r-base_r-pheatmap:e391eff8c37c62a5'}"

    input:
    path count_matrix_rds
    path samplesheet_csv
    val min_gs
    val min_mm
    val top_n
    val cor_threshold
    val pval_threshold
    val tomtype
    val networktype
    val deepsplit
    val reassignthreshold
    val mergecutheight
    val minmodulesize
    val sft_r2_threshold

    output:
    path "WGCNA", emit: wgcna_results
    path "WGCNA/vsd_object.rds", emit: vsd_object, optional: true
    path "WGCNA/unfiltered_counts.rds", emit: unfiltered_log2counts, optional: true
    path "WGCNA/filtered_counts.rds", emit: filtered_log2counts, optional: true
    path "WGCNA/normalized_counts.rds", emit: normalized_log2counts, optional: true
    path "WGCNA/all_traits_modules_hubgenes.rds", emit: modules_hubgenes
    path "WGCNA/soft_threshold_selection.pdf", emit: soft_threshold_plot, optional: true
    path "WGCNA/module_dendrogram.pdf", emit: module_dendrogram, optional: true
    path "WGCNA/module_trait_correlation_heatmap.pdf", emit: module_trait_heatmap, optional: true
    path "WGCNA/module_eigengenes_heatmap.pdf", emit: module_eigengenes_heatmap, optional: true
    path "WGCNA/module_analysis/**", emit: module_analysis_results, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def maxblocksize = task.memory.giga * 3500
    """
    run_wgcna_analysis.R ${count_matrix_rds} ${samplesheet_csv} ${min_gs} \\
    ${min_mm} ${top_n} ${cor_threshold} ${pval_threshold} ${maxblocksize} \\
    ${tomtype} ${networktype} ${deepsplit} ${reassignthreshold} \\
    ${mergecutheight} ${minmodulesize} ${sft_r2_threshold}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-wgcna: \$(Rscript -e "library(WGCNA); cat(as.character(packageVersion('WGCNA')))")
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
        r-pheatmap: \$(Rscript -e "library(pheatmap); cat(as.character(packageVersion('pheatmap')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p WGCNA
    mkdir -p WGCNA/module_analysis/sample_trait/csv
    mkdir -p WGCNA/module_analysis/sample_trait/plots

    touch WGCNA/unfiltered_counts.rds
    touch WGCNA/filtered_counts.rds
    touch WGCNA/normalized_counts.rds
    touch WGCNA/vsd_object.rds
    touch WGCNA/all_traits_modules_hubgenes.rds

    touch WGCNA/soft_threshold_selection.pdf
    touch WGCNA/module_dendrogram.pdf
    touch WGCNA/module_trait_correlation_heatmap.pdf
    touch WGCNA/module_eigengenes_heatmap.pdf

    touch WGCNA/module_analysis/sample_trait/csv/blue_genes.csv
    touch WGCNA/module_analysis/sample_trait/csv/blue_hub_genes.csv
    touch WGCNA/module_analysis/sample_trait/csv/blue_top${top_n}_hub_genes.csv
    touch WGCNA/module_analysis/sample_trait/plots/MM_vs_GS_blue.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-wgcna: \$(Rscript -e "library(WGCNA); cat(as.character(packageVersion('WGCNA')))")
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
        r-pheatmap: \$(Rscript -e "library(pheatmap); cat(as.character(packageVersion('pheatmap')))")
    END_VERSIONS
    """
}
