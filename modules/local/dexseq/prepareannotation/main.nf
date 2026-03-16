process DEXSEQ_PREPAREANNOTATION {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/06/068ebe7b3287bb1487cfd2276ce8bbb09be5921e7c4e3e5881813bd6d2792077/data'
        : 'community.wave.seqera.io/library/htseq:2.0.9--4a65a9021e1142a5'}"

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
