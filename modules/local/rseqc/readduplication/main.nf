process RSEQC_READDUPLICATION {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-c80ae8d0fe5685926c9bc673e400ff09a71844fd:29c8e89bc12d33b39e760c5ca3b1cfa087927580-0'
        : 'biocontainers/mulled-v2-c80ae8d0fe5685926c9bc673e400ff09a71844fd:e01414b01cd5729e641b84adaf1ce4fd6181bcb8-2'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.DupRate_plot.pdf"), emit: pdf
    tuple val(meta), path("*.DupRate_plot.r"), emit: rscript
    tuple val(meta), path("*.seq.DupRate.xls"), emit: seq_duplicaterate
    tuple val(meta), path("*.pos.DupRate.xls"), emit: pos_duplicaterate
    tuple val("${task.process}"), val('rseqc'), eval("read_duplication.py --version | sed 's/.*.py //'"), topic: versions, emit: versions_rseqc

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    read_duplication.py \\
        -i ${bam} \\
        ${args} \\
        -o ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.DupRate_plot.pdf
    touch ${prefix}.DupRate_plot.r
    touch ${prefix}.seq.DupRate.xls
    touch ${prefix}.pos.DupRate.xls
    """
}
