process RSEQC_JUNCTIONANNOTATION {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-c80ae8d0fe5685926c9bc673e400ff09a71844fd:29c8e89bc12d33b39e760c5ca3b1cfa087927580-0'
        : 'biocontainers/mulled-v2-c80ae8d0fe5685926c9bc673e400ff09a71844fd:e01414b01cd5729e641b84adaf1ce4fd6181bcb8-2'}"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(bed)

    output:
    tuple val(meta), path("*.xls"), emit: xls
    tuple val(meta), path("*.r"), emit: rscript
    tuple val(meta), path("*.log"), emit: log
    tuple val(meta), path("*.junction.bed"), optional: true, emit: junction_bed
    tuple val(meta), path("*.Interact.bed"), optional: true, emit: interact_bed
    tuple val(meta), path("*junction.pdf"), optional: true, emit: junction_pdf
    tuple val(meta), path("*events.pdf"), optional: true, emit: events_pdf
    tuple val("${task.process}"), val('rseqc'), eval("junction_annotation.py --version | sed 's/.*.py //'"), topic: versions, emit: versions_rseqc

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    junction_annotation.py \\
        -i ${bam} \\
        ${args} \\
        -o ${prefix} \\
        -r ${bed} \\
        2>| >(tee ${prefix}.junction_annotation.log >&2)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.junction.xls
    touch ${prefix}.junction_plot.r
    touch ${prefix}.junction.bed
    touch ${prefix}.Interact.bed
    touch ${prefix}.junction.pdf
    touch ${prefix}.events.pdf
    touch ${prefix}.junction_annotation.log
    """
}
