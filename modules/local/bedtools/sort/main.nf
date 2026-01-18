#!/usr/bin/env nextflow

process BEDTOOLS_SORT {
    tag "${bed.baseName}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--h13024bc_3'
        : 'biocontainers/bedtools:2.31.1--h13024bc_3'}"

    input:
    path bed

    output:
    path "*.sorted.bed", emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${bed.baseName}"
    """
    bedtools \\
        sort \\
        -i ${bed} \\
        ${args} \\
        > ${prefix}.sorted.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$( bedtools --version | sed 's/.*tools v//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${bed.baseName}"
    """
    touch ${prefix}.sorted.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$( bedtools --version | sed 's/.*tools v//' )
    END_VERSIONS
    """
}
