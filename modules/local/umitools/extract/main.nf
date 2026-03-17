process UMITOOLS_EXTRACT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/feb4b2af2e4e05b8a5e2f008e25b09508fff8fda2dd802c084617c6e16121210/data'
        : 'community.wave.seqera.io/library/umi_tools:1.1.4--8e2e1269f867d804'}"

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
