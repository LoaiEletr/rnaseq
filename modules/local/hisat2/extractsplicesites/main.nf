#!/usr/bin/env nextflow

process HISAT2_EXTRACTSPLICESITES {
    tag "${gtf.baseName}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/hisat2:2.2.1--h503566f_8'
        : 'biocontainers/hisat2:2.2.1--h503566f_8'}"

    input:
    path gtf

    output:
    path "*.splice_sites.tsv", emit: splice_sites
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${gtf.baseName}"
    """
    hisat2_extract_splice_sites.py \\
        ${gtf} \\
        ${args} \\
        > ${prefix}.splice_sites.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$( hisat2 --version | head -n 1 | awk '{print \$3}' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${gtf.baseName}"
    """
    touch ${prefix}.splice_sites.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$( hisat2 --version | head -n 1 | awk '{print \$3}' )
    END_VERSIONS
    """
}
