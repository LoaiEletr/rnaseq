#!/usr/bin/env nextflow

process FASTQC {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0'
        : 'biocontainers/fastqc:0.12.1--hdfd78af_0'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"), emit: zip
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fastqc \\
        ${args} \\
        --threads ${task.cpus} \\
        ${reads} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed 's/.*v//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def touch_fastqc_output = meta.single_end ? "touch ${prefix}_fastqc.html ; touch ${prefix}_fastqc.zip" : "touch ${prefix}_1_fastqc.html  ; touch ${prefix}_2_fastqc.html ; touch ${prefix}_1_fastqc.zip  ; touch ${prefix}_2_fastqc.zip"
    """
    ${touch_fastqc_output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed 's/.*v//' )
    END_VERSIONS
    """
}
