#!/usr/bin/env nextflow

process TX2GENE {
    tag "${species_name}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f5/f58661c83be1d4605d0902036012cf776c0744f8fdc1ec3c5d8dcaed038b0553/data'
        : 'community.wave.seqera.io/library/bioconductor-biomart_r-base:f9e5215e8454b7b6'}"

    input:
    val species_name

    output:
    path "*.tsv", emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    generate_tx2gene.R ${species_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-biomart: \$(Rscript -e "library(biomaRt); cat(as.character(packageVersion('biomaRt')))")
    END_VERSIONS
    """

    stub:
    """
    touch tx2gene.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-biomart: \$(Rscript -e "library(biomaRt); cat(as.character(packageVersion('biomaRt')))")
    END_VERSIONS
    """
}
