process GXF2BED {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/93/93bc201ef43d1be7fd3cab3635f7bf7dd7bcc55ad753c726bf3b7bcbb3f6f802/data'
        : 'community.wave.seqera.io/library/gxf2bed:0.2.7--1aa890a939753630'}"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val("${task.process}"), val('gxf2bed'), eval("gxf2bed --version | grep 'gxf2bed [0-9]' | sed 's/.*ed //'"), topic: versions, emit: versions_gxf2bed

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gxf2bed \\
        ${args} \\
        -i ${gtf} \\
        -o ${prefix}.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    """
}
