#!/usr/bin/env nextflow

process UMITOOLS_DEDUP {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/umi_tools:1.1.4--py39hf95cd2a_2'
        : 'biocontainers/umi_tools:1.1.4--py39hf95cd2a_2'}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.log"), emit: log
    tuple val(meta), path("*edit_distance.tsv"), emit: tsv_edit_distance
    tuple val(meta), path("*per_umi.tsv"), emit: tsv_per_umi
    tuple val(meta), path("*per_position.tsv"), emit: tsv_umi_per_position
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bam_in = meta.single_end ? "" : "--paired"
    """
    umi_tools \\
        dedup \\
        ${args} \\
        -I ${bam} \\
        -L ${prefix}.umi_dedup.log \\
        -S ${prefix}.umi_dedup.bam \\
        ${bam_in} \\
        --extract-umi-method=read_id \\
        --umi-separator=":" \\
        --output-stats=${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        umi_tools: \$( umi_tools --version | sed 's/.*version: //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.umi_dedup.bam
    touch ${prefix}.umi_dedup.log
    touch ${prefix}.edit_distance.tsv
    touch ${prefix}.per_umi.tsv
    touch ${prefix}.per_position.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        umi_tools: \$( umi_tools --version | sed 's/.*version: //' )
    END_VERSIONS
    """
}
