process UMITOOLS_EXTRACT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/umi_tools:1.1.6--py39hbcbf7aa_0'
        : 'biocontainers/umi_tools:1.1.6--py310h1fe012e_0'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('umitools'), eval("umi_tools --version | sed 's/.*version: //'"), topic: versions, emit: versions_umitools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fastq_in = meta.single_end ? "-I ${reads}" : "-I ${reads[0]} --read2-in ${reads[1]}"
    def fastq_out = meta.single_end ? "-S ${prefix}.umi_extract.fastq.gz" : "-S ${prefix}.umi_extract_1.fastq.gz --read2-out ${prefix}.umi_extract_2.fastq.gz"
    """
    umi_tools extract \\
        ${args} \\
        ${fastq_in} \\
        ${fastq_out} \\
        --log=${prefix}.umi_extract.log
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gzip_fastq_out = meta.single_end ? "echo '' | gzip > ${prefix}.umi_extract.fastq.gz" : "echo '' | gzip > ${prefix}.umi_extract_1.fastq.gz ; echo '' | gzip > ${prefix}.umi_extract_2.fastq.gz"
    """
    ${gzip_fastq_out}
    touch ${prefix}.umi_extract.log
    """
}
