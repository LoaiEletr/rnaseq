#!/usr/bin/env nextflow

process LIMMA_DE {
    tag "${voom_object}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd1284a73c2a23f3cef38ec39f1695005012adab9c3ecc4ba56ed33a0b4daa14/data'
        : 'community.wave.seqera.io/library/bioconductor-limma_r-base:2579913b34b8b4da'}"

    input:
    path voom_object
    path samplesheet
    val pvalue_threshold
    val logfc_threshold
    val top_genes

    output:
    path "ebfit.rds", emit: ebfit_rds
    path "de_expression_values.rds", emit: de_expression_rds
    path "de_genes_full_table.rds", emit: de_genes_rds
    path "unfiltered_deg_results.rds", emit: unfiltered_deg_rds
    path "top_${top_genes}_genes.rds", emit: topgenes_rds
    path "top_${top_genes}_genes.csv", emit: topgenes_csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_limma_de.R ${voom_object} ${samplesheet} ${pvalue_threshold} ${logfc_threshold} ${top_genes}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-limma: \$(Rscript -e "library(limma); cat(as.character(packageVersion('limma')))" 2>/dev/null || echo "NA")
    END_VERSIONS
    """

    stub:
    """
    touch ebfit.rds
    touch de_expression_values.rds
    touch de_genes_full_table.rds
    touch unfiltered_deg_results.rds
    touch top_${top_genes}_genes.rds
    touch top_${top_genes}_genes.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-limma: \$(Rscript -e "library(limma); cat(as.character(packageVersion('limma')))" 2>/dev/null || echo "NA")
    END_VERSIONS
    """
}
