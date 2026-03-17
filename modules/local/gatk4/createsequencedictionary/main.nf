process GATK4_CREATESEQUENCEDICTIONARY {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.dict"), emit: dict
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def avail_mem = task.memory.mega
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        CreateSequenceDictionary \\
        ${args} \\
        --REFERENCE ${fasta} \\
        --URI ${fasta} \\
        --TMP_DIR .
    """

    stub:
    """
    touch ${fasta.baseName}.dict
    """
}
