#!/usr/bin/env nextflow

process MULTIQC {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/multiqc:1.26--pyhdfd78af_0'
        : 'biocontainers/multiqc:1.33--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(multiqc_files)
    tuple val(meta2), path(multiqc_config)

    output:
    tuple val(meta), path("*.html"), emit: report
    tuple val(meta), path("*_data"), emit: data
    tuple val(meta), path("*_plots"), optional: true, emit: plots
    tuple val("${task.process}"), val('multiqc'), eval("multiqc --version | sed 's/.*version //'"), emit: versions_multiqc

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def config = multiqc_config ? "--config ${multiqc_config}" : ''
    """
    multiqc \\
        --force \\
        ${args} \\
        ${config} \\
        ${multiqc_files} \\
    """

    stub:
    """
    mkdir multiqc_data
    mkdir multiqc_plots
    touch multiqc_report.html
    """
}
