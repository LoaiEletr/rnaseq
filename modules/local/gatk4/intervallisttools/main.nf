process GATK4_INTERVALLISTTOOLS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(intervals)

    output:
    tuple val(meta), path("*_split/*/*.interval_list"), emit: interval_list
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = task.memory.mega
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        IntervalListTools \\
        ${args} \\
        --INPUT ${intervals} \\
        --OUTPUT ${prefix}_split \\
        --TMP_DIR .

    COUNT=1
    # The glob corresponds to: folder_split/subfolder/filename.interval_list
    for interval in *_split/*/*.interval_list; do
        DIR=\$(dirname "\${interval}")
        BASE=\$(basename "\${interval}")
        mv "\${interval}" "\${DIR}/\${COUNT}\${BASE}"
        COUNT=\$((COUNT + 1))
    done
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_split/temp_0001_of_4
    mkdir -p ${prefix}_split/temp_0002_of_4
    mkdir -p ${prefix}_split/temp_0003_of_4
    mkdir -p ${prefix}_split/temp_0004_of_4
    touch ${prefix}_split/temp_0001_of_4/1scattered.interval_list
    touch ${prefix}_split/temp_0002_of_4/2scattered.interval_list
    touch ${prefix}_split/temp_0003_of_4/3scattered.interval_list
    touch ${prefix}_split/temp_0004_of_4/4scattered.interval_list
    """
}
