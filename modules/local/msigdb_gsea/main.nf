#!/usr/bin/env nextflow

process MSIGDB_GSEA {
    tag "${unfiltered_deg}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9e/9e4333b647c7183896fe4cd6054c7733a9a5b93771d50acda8a5169436591e02/data'
        : 'community.wave.seqera.io/library/bioconductor-clusterprofiler_bioconductor-enrichplot_r-base_r-dplyr_pruned:26af7e9ffc04c073'}"

    input:
    path unfiltered_deg
    val pvalue_threshold
    val species_name
    val msigdb_categories
    val nes_threshold
    val padj_gsea
    val ntop_processes
    val rank_method

    output:
    path "GSEA_results", emit: gsea_results, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_msigdb_gsea.R ${unfiltered_deg} ${pvalue_threshold} ${species_name} ${msigdb_categories} ${nes_threshold} ${padj_gsea} ${ntop_processes} ${rank_method}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-clusterprofiler: \$(Rscript -e "library(clusterProfiler); cat(as.character(packageVersion('clusterProfiler')))")
        bioconductor-enrichplot: \$(Rscript -e "library(enrichplot); cat(as.character(packageVersion('enrichplot')))")
        r-ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        r-dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        r-tibble: \$(Rscript -e "library(tibble); cat(as.character(packageVersion('tibble')))")
        r-msigdbr: \$(Rscript -e "library(msigdbr); cat(as.character(packageVersion('msigdbr')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir GSEA_results

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-clusterprofiler: \$(Rscript -e "library(clusterProfiler); cat(as.character(packageVersion('clusterProfiler')))")
        bioconductor-enrichplot: \$(Rscript -e "library(enrichplot); cat(as.character(packageVersion('enrichplot')))")
        r-ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        r-dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        r-tibble: \$(Rscript -e "library(tibble); cat(as.character(packageVersion('tibble')))")
        r-msigdbr: \$(Rscript -e "library(msigdbr); cat(as.character(packageVersion('msigdbr')))")
    END_VERSIONS
    """
}
