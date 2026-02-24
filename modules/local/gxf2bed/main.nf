process GXF2BED {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gxf2bed:0.2.4--ha6fb395_0'
        : 'biocontainers/gxf2bed:0.2.7--ha6fb395_0'}"

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
