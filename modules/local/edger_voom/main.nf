#!/usr/bin/env nextflow

process EDGER_VOOM {
    tag "${count_matrix_rds}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/60/60b30077e86da5ad7048eb82ef265e1c525dbc3ca9d284f5c1bac530cdf33549/data'
        : 'community.wave.seqera.io/library/bioconductor-edger_r-base:4.4.0--11c47d95799573a7'}"

    input:
    path count_matrix_rds
    path samplesheet_csv

    output:
    path "unfiltered_log2cpm.rds", emit: unfiltered_log2cpm
    path "filtered_log2cpm.rds", emit: filtered_log2cpm
    path "normalized_log2cpm.rds", emit: normalized_log2cpm
    path "voom_transformed.rds", emit: voom_transformed
    path "voom_plot.pdf", emit: voom_plot
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_edger_voom.R ${count_matrix_rds} ${samplesheet_csv}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-edger: \$(Rscript -e "library(edgeR); cat(as.character(packageVersion('edgeR')))")
    END_VERSIONS
    """

    stub:
    """
    touch unfiltered_log2cpm.rds
    touch filtered_log2cpm.rds
    touch normalized_log2cpm.rds
    touch voom_transformed.rds
    touch voom_plot.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-edger: \$(Rscript -e "library(edgeR); cat(as.character(packageVersion('edgeR')))")
    END_VERSIONS
    """
}
