process SEQKIT_STATS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'
        : 'biocontainers/seqkit:2.12.0--he881be0_1'}"

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
