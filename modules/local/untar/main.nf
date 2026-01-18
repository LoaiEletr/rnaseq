#!/usr/bin/env nextflow

process UNTAR {
    tag "${archive}"
    label 'process_low'

    input:
    path archive

    output:
    path ("${prefix}"), emit: untar
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = archive.toString() - '.tar.gz'
    """
    tar \\
        -xzvf \\
        ${args} \\
        ${archive} \\
        --one-top-level=${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        untar: \$( tar --version | head -n 1 | sed 's/.*) //' )
    END_VERSIONS
    """

    stub:
    prefix = archive.toString() - '.tar.gz'
    """
    touch ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        untar: \$( tar --version | head -n 1 | sed 's/.*) //' )
    END_VERSIONS
    """
}
