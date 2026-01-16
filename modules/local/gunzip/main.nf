#!/usr/bin/env nextflow

process GUNZIP {
    tag "${archive}"
    label 'process_low'

    input:
    tuple val(meta), path(archive)

    output:
    tuple val(meta), path("${gunzip}"), emit: gunzip
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    gunzip = archive.toString() - '.gz'
    """
    gzip \\
        -cd \\
        ${args} \\
        -c ${archive} \\
        > ${gunzip}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunzip: \$( gunzip --version | grep "gunzip" | awk '{print \$3}' )
    END_VERSIONS
    """

    stub:
    gunzip = archive.toString() - '.gz'
    """
    touch ${gunzip}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunzip: \$( gunzip --version | grep "gunzip" | awk '{print \$3}' )
    END_VERSIONS
    """
}
