process SEQKIT_STATS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/85/85b40b925e4d4a62f9b833bbb0646d7ea6cf53d8a875e3055f90da757d7ccd27/data'
        : 'community.wave.seqera.io/library/seqkit:2.12.0--5e60eb93d3a59212'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), env("AVG_READLENGTH"), emit: avg_length
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/seqkit v//'"), topic: versions, emit: versions_seqkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def fastq_in = meta.single_end ? "${reads}" : "${reads[0]}"
    """
    AVG_READLENGTH=\$(seqkit stats ${fastq_in} | awk 'NR==2 {printf "%d\\n", \$8+0.5}')
    """

    stub:
    """
    AVG_READLENGTH=100
    """
}
