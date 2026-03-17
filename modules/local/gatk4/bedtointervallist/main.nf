process GATK4_BEDTOINTERVALLIST {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(bed)
    tuple val(meta2), path(dict)

    output:
    tuple val(meta), path("*.interval_list"), emit: interval_list
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = task.memory.mega
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        BedToIntervalList \\
        ${args} \\
        --INPUT ${bed} \\
        --OUTPUT ${prefix}.interval_list \\
        --SEQUENCE_DICTIONARY ${dict} \\
        --TMP_DIR .
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.interval_list
    """
}
