#!/usr/bin/env nextflow

process PICARD_MARKDUPLICATES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/08/0861295baa7c01fc593a9da94e82b44a729dcaf8da92be8e565da109aa549b25/data'
        : 'community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.metrics.txt"), emit: metrics
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = task.memory.mega
    """
    picard \\
        -Xmx${avail_mem}M \\
        MarkDuplicates \\
        ${args} \\
        --INPUT ${bam} \\
        --OUTPUT ${prefix}.bam \\
        --METRICS_FILE ${prefix}.MarkDuplicates.metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$( picard MarkDuplicates --version 2>&1 | sed -n 's/^Version:*//p' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.MarkDuplicates.metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$( picard MarkDuplicates --version 2>&1 | sed -n 's/^Version:*//p' )
    END_VERSIONS
    """
}
