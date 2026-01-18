#!/usr/bin/env nextflow

process CUTADAPT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/cutadapt:5.0--py39hbcbf7aa_0'
        : 'biocontainers/cutadapt:5.2--py310h1fe012e_0'}"

    input:
    tuple val(meta), path(reads)
    path adapter

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val(meta), path("*.json"), emit: json
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def adapter_in = meta.single_end ? "-a file:${adapter}" : "-a file:${adapter} -A file:${adapter}"
    def fastq_in = meta.single_end ? "-o ${prefix}.cutadapt.fastq.gz" : "-o ${prefix}.cutadapt_1.fastq.gz -p ${prefix}.cutadapt_2.fastq.gz"
    def fastq_out = meta.single_end ? "${reads}" : "${reads[0]} ${reads[1]}"
    """
    cutadapt \\
        ${args} \\
        ${args2} \\
        --json=${prefix}.cutadapt.json \\
        ${adapter_in} \\
        ${fastq_out} \\
        ${fastq_in}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$( cutadapt --version )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gzip_fastq_out = meta.single_end ? "echo '' | gzip > ${prefix}.cutadapt.fastq.gz" : "echo '' | gzip > ${prefix}.cutadapt_1.fastq.gz ; echo '' | gzip > ${prefix}.cutadapt_2.fastq.gz"
    """
    ${gzip_fastq_out}
    touch ${prefix}.cutadapt.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$( cutadapt --version )
    END_VERSIONS
    """
}
