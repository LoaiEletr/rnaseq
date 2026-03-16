process UMITOOLS_DEDUP {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/38/38136a834cc8aecb5a477640910df9472eb83e90118059bd1a0459aec685b731/data'
        : 'community.wave.seqera.io/library/umi_tools:1.1.6--2c79f13606642e4c'}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.log"), emit: log
    tuple val(meta), path("*edit_distance.tsv"), emit: tsv_edit_distance
    tuple val(meta), path("*per_umi.tsv"), emit: tsv_per_umi
    tuple val(meta), path("*per_position.tsv"), emit: tsv_umi_per_position
    tuple val("${task.process}"), val('umitools'), eval("umi_tools --version | sed 's/.*version: //'"), topic: versions, emit: versions_umitools

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
        --output-stats=${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.umi_dedup.bam
    touch ${prefix}.umi_dedup.log
    touch ${prefix}.edit_distance.tsv
    touch ${prefix}.per_umi.tsv
    touch ${prefix}.per_position.tsv
    """
}
