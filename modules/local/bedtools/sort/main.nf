process BEDTOOLS_SORT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7d/7df273d12f0c4d8539440b68876edf39b739cb78bb806418c5b5d057fe11bdbd/data'
        : 'community.wave.seqera.io/library/bedtools:2.31.1--7c4ce4cb07c09ee4'}"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.sorted.bed"), emit: bed
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed 's/.*tools v//'"), topic: versions, emit: versions_bedtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bedtools \\
        sort \\
        -i ${bed} \\
        ${args} \\
        > ${prefix}.sorted.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sorted.bed
    """
}
