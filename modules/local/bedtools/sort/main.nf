process BEDTOOLS_SORT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--h13024bc_3'
        : 'biocontainers/bedtools:2.31.1--h13024bc_3'}"

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
