#!/usr/bin/env nextflow

process DEXSEQ_PREPAREANNOTATION {
    tag "${gtf}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/htseq:2.0.5--py39h8931b72_3'
        : 'biocontainers/htseq:2.0.9--py312h8f4af18_0'}"

    input:
    path gtf
    val aggregation

    output:
    path "*.gff", emit: gff
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "DEXSeq"

    def aggregation_flag = aggregation ? '' : '-r no'
    """
    dexseq_prepare_annotation.py \\
        ${gtf} \\
        ${prefix}.gff \\
        ${args} \\
        ${aggregation_flag}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htseq: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('htseq').version)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "DEXSeq"
    """
    touch ${prefix}.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htseq: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('htseq').version)")
    END_VERSIONS
    """
}
