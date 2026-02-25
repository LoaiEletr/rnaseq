process SUBREAD_FEATURECOUNTS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/subread:2.0.8--h577a1d6_0'
        : 'biocontainers/subread:2.1.1--h577a1d6_0'}"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("*.counts.txt"), emit: counts
    tuple val(meta), path("*.summary"), emit: summary
    tuple val("${task.process}"), val('subread'), eval("featureCounts -v 2>&1 | sed 's/featureCounts v//'"), topic: versions, emit: versions_subread

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bam_in = meta.single_end ? "${bam}" : "-p --countReadPairs ${bam}"
    strand_flag = ""
    if (meta.lib_type in ["IU", "U", "MU", "OU"]) {
        strand_flag = "-s 0"
    }
    else if (meta.lib_type in ["SR", "ISR", "MSR", "OSR", "reverse"]) {
        strand_flag = "-s 2"
    }
    else if (meta.lib_type in ["SF", "ISF", "MSF", "OSF", "forward"]) {
        strand_flag = "-s 1"
    }
    else {
        error("Unknown strandness type: ${meta.lib_type}")
    }
    """
    featureCounts \\
        ${args} \\
        ${strand_flag} \\
        -t exon \\
        -g gene_id \\
        -a ${gtf} \\
        -o ${prefix}.counts.txt \\
        ${bam_in}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.counts.txt
    touch ${prefix}.counts.txt.summary
    """
}
