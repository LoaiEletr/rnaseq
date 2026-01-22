#!/usr/bin/env nextflow

process TXIMPORT {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b8/b8617b62baac3163f0458b3f72ede655a6376c5b713f7d7b7507c6dbe7209f9e/data'
        : 'community.wave.seqera.io/library/bioconductor-rhdf5_bioconductor-tximport_r-base_r-jsonlite:0454e92595745a41'}"

    input:
    path quant_dirs
    val quant_type
    path tx2gene_table

    output:
    path "*.rds", emit: rds
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_tximport.R ${quant_dirs} ${quant_type} ${tx2gene_table}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-tximport: \$(Rscript -e "library(tximport); cat(as.character(packageVersion('tximport')))")
        bioconductor-rhdf5: \$(Rscript -e "library(rhdf5); cat(as.character(packageVersion('rhdf5')))")
    END_VERSIONS
    """

    stub:
    """
    touch tximport_results.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-tximport: \$(Rscript -e "library(tximport); cat(as.character(packageVersion('tximport')))")
        bioconductor-rhdf5: \$(Rscript -e "library(rhdf5); cat(as.character(packageVersion('rhdf5')))")
    END_VERSIONS
    """
}
