#!/usr/bin/env nextflow

process MASIGPRO_TIMECOURSE {
    tag "${count_matrix_rds}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a1/a15ce4db557dad1855f663faeb430ad604b0cf15e64b25f63663199d5bd5cc01/data'
        : 'community.wave.seqera.io/library/bioconductor-edger_bioconductor-masigpro_r-base_r-dplyr_pruned:6895f6766e88b3d8'}"

    input:
    path count_matrix_rds
    path samplesheet_csv
    val pvalue_threshold
    val rsq_threshold
    val cluster_method

    output:
    path "masigpro_results", emit: masigpro_results
    path "masigpro_results/unfiltered_log2cpm.rds", emit: unfiltered_log2cpm
    path "masigpro_results/filtered_log2cpm.rds", emit: filtered_log2cpm
    path "masigpro_results/normalized_log2cpm.rds", emit: normalized_log2cpm
    path "masigpro_results/*clusterplots.pdf", emit: cluster_plots
    path "masigpro_results/*median_expression.csv", emit: cluster_median_expression
    path "masigpro_results/*statistics.csv", emit: summary_statistics
    path "masigpro_results/*assignments.csv", emit: cluster_assignments
    path "masigpro_results/expression_matrix_siggenes.rds", emit: significant_expression
    path "masigpro_results/clusters_list.rds", emit: clusters_list
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_masigpro.R ${count_matrix_rds} ${samplesheet_csv} ${pvalue_threshold} ${rsq_threshold} ${cluster_method}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-masigpro: \$(Rscript -e "library(maSigPro); cat(as.character(packageVersion('maSigPro')))")
        r-dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        r-tibble: \$(Rscript -e "library(tibble); cat(as.character(packageVersion('tibble')))")
        bioconductor-edger: \$(Rscript -e "library(edgeR); cat(as.character(packageVersion('edgeR')))")
        r-ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        r-reshape2: \$(Rscript -e "library(reshape2); cat(as.character(packageVersion('reshape2')))")
        r-patchwork: \$(Rscript -e "library(patchwork); cat(as.character(packageVersion('patchwork')))")
        r-mclust: \$(Rscript -e "library(mclust); cat(as.character(packageVersion('mclust')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p masigpro_results
    touch masigpro_results/cluster_assignments.csv
    touch masigpro_results/cluster_median_expression.csv
    touch masigpro_results/summary_statistics.csv
    touch masigpro_results/expression_matrix_siggenes.rds
    touch masigpro_results/clusters_list.rds
    touch masigpro_results/unfiltered_log2cpm.rds
    touch masigpro_results/filtered_log2cpm.rds
    touch masigpro_results/normalized_log2cpm.rds
    touch masigpro_results/masigpro_clusterplots.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-masigpro: \$(Rscript -e "library(maSigPro); cat(as.character(packageVersion('maSigPro')))")
        r-dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        r-tibble: \$(Rscript -e "library(tibble); cat(as.character(packageVersion('tibble')))")
        bioconductor-edger: \$(Rscript -e "library(edgeR); cat(as.character(packageVersion('edgeR')))")
        r-ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        r-reshape2: \$(Rscript -e "library(reshape2); cat(as.character(packageVersion('reshape2')))")
        r-patchwork: \$(Rscript -e "library(patchwork); cat(as.character(packageVersion('patchwork')))")
        r-mclust: \$(Rscript -e "library(mclust); cat(as.character(packageVersion('mclust')))")
    END_VERSIONS
    """
}
