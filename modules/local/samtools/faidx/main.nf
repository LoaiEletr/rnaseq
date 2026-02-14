process SAMTOOLS_FAIDX {
    tag "${fasta}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.21--h96c455f_1'
        : 'biocontainers/samtools:1.23--h96c455f_0'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fai"), emit: fai
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | head -n 1 | sed 's/samtools //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools \\
        faidx \\
        ${args} \\
        ${fasta}
    """

    stub:
    """
    touch ${fasta}.fai
    """
}
