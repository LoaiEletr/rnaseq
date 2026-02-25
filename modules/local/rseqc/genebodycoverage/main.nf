process RSEQC_GENEBODYCOVERAGE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-c80ae8d0fe5685926c9bc673e400ff09a71844fd:29c8e89bc12d33b39e760c5ca3b1cfa087927580-0'
        : 'biocontainers/mulled-v2-c80ae8d0fe5685926c9bc673e400ff09a71844fd:e01414b01cd5729e641b84adaf1ce4fd6181bcb8-2'}"

    input:
    tuple val(meta), path(bam_files)
    tuple val(meta2), path(bai_files)
    tuple val(meta3), path(bed)

    output:
    tuple val(meta), path("*.geneBodyCoverage.r"), emit: rscript
    tuple val(meta), path("*.geneBodyCoverage.curves.pdf"), emit: curves
    tuple val(meta), path("*.geneBodyCoverage.heatMap.pdf"), optional: true, emit: heatmap
    tuple val(meta), path("*.geneBodyCoverage.txt"), emit: txt
    tuple val(meta), path("log.txt"), emit: log
    tuple val("${task.process}"), val('rseqc'), eval("geneBody_coverage.py --version | sed 's/.*.py //'"), topic: versions, emit: versions_rseqc

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "genebodycoverage"
    """
    geneBody_coverage.py \\
        -r ${bed} \\
        -i ${bam_files.join(',')} \\
        -o ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "genebodycoverage"
    """
    touch ${prefix}.geneBodyCoverage.r
    touch ${prefix}.geneBodyCoverage.curves.pdf
    touch ${prefix}.geneBodyCoverage.heatMap.pdf
    touch ${prefix}.geneBodyCoverage.txt
    touch log.txt
    """
}
