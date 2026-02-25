process DEXSEQ_PREPAREANNOTATION {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/htseq:2.0.5--py39h8931b72_3'
        : 'biocontainers/htseq:2.0.9--py312h8f4af18_0'}"

    input:
    tuple val(meta), path(gtf)
    val aggregation

    output:
    tuple val(meta), path("*.gff"), emit: gff
    tuple val("${task.process}"), val('htseq'), eval("python -c 'from importlib.metadata import version; print(version(\"HTSeq\"))'"), topic: versions, emit: versions_htseq

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "DEXSeq"

    def aggregation_flag = aggregation ? '' : '-r no'
    """
    dexseq_prepare_annotation.py \\
        ${gtf} \\
        ${prefix}.gff \\
        ${args} \\
        ${aggregation_flag}
    """

    stub:
    def prefix = task.ext.prefix ?: "DEXSeq"
    """
    touch ${prefix}.gff
    """
}
