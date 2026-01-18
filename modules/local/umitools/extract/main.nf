#!/usr/bin/env nextflow

process UMITOOLS_EXTRACT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/umi_tools:1.1.6--py39hbcbf7aa_0'
        : 'biocontainers/umi_tools:1.1.6--py310h1fe012e_0'}"

    input:
    tuple val(meta), path(reads)
    val lib_kit

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def kit_map = [
        "quantseq": [
            umitools_extract_method: "regex",
            umitools_bc_pattern: "'^(?P<umi_1>.{6})(?P<discard_1>.{4}).*'",
        ],
        "corall": [
            umitools_extract_method: "regex",
            umitools_bc_pattern: "'^(?P<umi_1>.{12}).*'",
        ],
        "takara": [
            umitools_extract_method: "regex",
            umitools_bc_pattern2: "'^(?P<umi_1>.{8})(?P<discard_1>.{6}).*'",
        ],
    ]

    def config = kit_map["${lib_kit}"]
    def extract_method = config.umitools_extract_method
    def bc_pattern = config.umitools_bc_pattern ? "--bc-pattern ${config.umitools_bc_pattern}" : ""
    def bc_pattern2 = config.umitools_bc_pattern2 ? "--bc-pattern2 ${config.umitools_bc_pattern2}" : ""
    def fastq_in = meta.single_end ? "-I ${reads}" : "-I ${reads[0]} --read2-in ${reads[1]}"
    def fastq_out = meta.single_end ? "-S ${prefix}.umi_extract.fastq.gz" : "-S ${prefix}.umi_extract_1.fastq.gz --read2-out ${prefix}.umi_extract_2.fastq.gz"
    """
    umi_tools extract \\
        --extract-method=${extract_method} \\
        --umi-separator=":" \\
        ${args} \\
        ${bc_pattern} \\
        ${bc_pattern2} \\
        ${fastq_in} \\
        ${fastq_out} \\
        --log=${prefix}.umi_extract.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        umi_tools: \$( umi_tools --version | sed 's/.*version: //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.umi_extract.fastq.gz
    touch ${prefix}.umi_extract.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        umi_tools: \$( umi_tools --version | sed 's/.*version: //' )
    END_VERSIONS
    """
}
