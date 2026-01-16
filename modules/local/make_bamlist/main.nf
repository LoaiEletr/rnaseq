#!/usr/bin/env nextflow

process MAKE_BAMLIST {
    tag "${meta.condition}"
    label 'process_low'

    input:
    tuple val(meta), path(bam_files)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.condition}"
    """
    echo \\
        ${bam_files.join(',')} \\
        > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | awk '{print \$4}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.condition}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | awk '{print \$4}')
    END_VERSIONS
    """
}
