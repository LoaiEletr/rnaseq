#!/usr/bin/env nextflow

process HISAT2_BUILD {
    tag "${ref.baseName}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/hisat2:2.2.1--h503566f_8'
        : 'biocontainers/hisat2:2.2.1--h503566f_8'}"

    input:
    path ref
    path splice_sites
    path exons_sites

    output:
    path "index", emit: index
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${ref.getSimpleName()}"
    """
    mkdir index
    hisat2-build \\
        --ss ${splice_sites} \\
        ${args} \\
        --exon ${exons_sites} \\
        ${ref} \\
        index/${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$( hisat2 --version | head -n 1 | awk '{print \$3}' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${ref.getSimpleName()}"
    """
    mkdir index
    touch "index/${prefix}."{1..8}.ht2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$( hisat2 --version | head -n 1 | awk '{print \$3}' )
    END_VERSIONS
    """
}
