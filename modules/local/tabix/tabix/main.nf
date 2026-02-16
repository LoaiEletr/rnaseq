process TABIX_TABIX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/02/02846e0365044645f2c1a35c0867431eec698cf3b0acd338df01703825627b28/data'
        : 'community.wave.seqera.io/library/htslib:1.23--7c40560d4c2ce27e'}"

    input:
    tuple val(meta), path(tab)

    output:
    tuple val(meta), path("*.tbi"), emit: index
    tuple val("${task.process}"), val('tabix'), eval("tabix -h 2>&1 | grep -oP 'Version:\\s*\\K[^\\s]+'"), topic: versions, emit: versions_tabix

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    tabix \\
        ${args} \\
        --threads ${task.cpus} \\
        ${tab}
    """

    stub:
    """
    touch ${tab}.tbi
    """
}
