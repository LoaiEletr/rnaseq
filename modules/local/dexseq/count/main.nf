process DEXSEQ_COUNT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/htseq:2.0.5--py39h8931b72_3'
        : 'biocontainers/htseq:2.0.9--py312h8f4af18_0'}"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(gff)
    val alignment_quality

    output:
    tuple val(meta), path("*.counts.txt"), emit: counts
    tuple val("${task.process}"), val('htseq'), eval("python -c 'from importlib.metadata import version; print(version(\"HTSeq\"))'"), topic: versions, emit: versions_htseq

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_type = meta.single_end ? '' : '-p yes'
    def alignment_quality_flag = "-a ${alignment_quality}"
    strand_flag = ""

    if (meta.lib_type in ["IU", "U", "MU", "OU"]) {
        strand_flag = "-s no"
    }
    else if (meta.lib_type in ["SR", "ISR", "MSR", "OSR", "reverse"]) {
        strand_flag = "-s reverse"
    }
    else if (meta.lib_type in ["SF", "ISF", "MSF", "OSF", "forward"]) {
        strand_flag = "-s yes"
    }
    else {
        error("Unknown strandness type: ${meta.lib_type}")
    }
    """
    dexseq_count.py \\
        ${gff} \\
        ${read_type} \\
        ${args} \\
        -f bam ${bam} \\
        -r pos \\
        ${prefix}.counts.txt \\
        ${alignment_quality_flag} \\
        ${strand_flag}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.counts.txt
    """
}
