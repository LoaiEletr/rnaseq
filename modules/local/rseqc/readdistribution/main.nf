process RSEQC_READDISTRIBUTION {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1d/1db6950626c14db8a5b5a80089c259774a693fc3a5946d1bf169d19b11f7bccb/data'
        : 'community.wave.seqera.io/library/rseqc_r-base:a2f5852c8ab06f43'}"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(bed)

    output:
    tuple val(meta), path("*.read_distribution.txt"), emit: txt
    tuple val("${task.process}"), val('rseqc'), eval("read_distribution.py --version | sed 's/.*.py //'"), topic: versions, emit: versions_rseqc

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    read_distribution.py \\
        -r ${bed} \\
        ${args} \\
        -i ${bam} \\
        > ${prefix}.read_distribution.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.read_distribution.txt
    """
}
