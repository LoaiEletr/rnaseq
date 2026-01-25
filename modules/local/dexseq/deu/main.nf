#!/usr/bin/env nextflow

process DEXSEQ_DEU {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b0/b0fba3f6ff907ea4d9acc1d639f02bec14fa6a7b771ec083e34d8a63456da2bd/data'
        : 'community.wave.seqera.io/library/bioconductor-dexseq_r-base:23680634f04cf138'}"

    input:
    path counts_txt
    path samplesheet_csv
    path flattened_gff
    val pvalue_threshold
    val logfc_threshold
    val top_n
    val min_exonlength

    output:
    path "DEXSeq_table.csv", optional: true, emit: results_table
    path "significant_DEUs.csv", optional: true, emit: significant_results
    path "DEXSeq_MA_plot.pdf", optional: true, emit: ma_plot
    path "top_DEU_plots", optional: true, emit: gene_plots
    path "significant_DEU_gene_ids.rds", optional: true, emit: significant_genes
    path "DEXSeqReport", optional: true, emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_dexseq_deu.R ${samplesheet_csv} ${flattened_gff} ${pvalue_threshold} ${logfc_threshold} ${top_n} ${min_exonlength}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-dexseq: \$(Rscript -e "library(DEXSeq); cat(as.character(packageVersion('DEXSeq')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir top_DEU_plots
    mkdir DEXSeqReport
    touch DEXSeq_table.csv
    touch significant_DEUs.csv
    touch DEXSeq_MA_plot.pdf
    touch significant_DEU_gene_ids.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-dexseq: \$(Rscript -e "library(DEXSeq); cat(as.character(packageVersion('DEXSeq')))")
    END_VERSIONS
    """
}
