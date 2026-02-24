process KALLISTO_INDEX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/kallisto:0.51.1--h2b92561_2'
        : 'biocontainers/kallisto:0.51.1--h2b92561_2'}"

    input:
    tuple val(meta), path(transcriptome)

    output:
    tuple val(meta), path("*.idx"), emit: index
    tuple val("${task.process}"), val('kallisto'), eval("kallisto version | sed 's/.*ion //'"), topic: versions, emit: versions_kallisto

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.argg ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    kallisto \\
        index \\
        -t ${task.cpus} \\
        ${args} \\
        -i ${prefix}.idx \\
        ${transcriptome}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.idx
    """
}
