#!/usr/bin/env nextflow

process RSEQC_INNERDISTANCE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-c80ae8d0fe5685926c9bc673e400ff09a71844fd:29c8e89bc12d33b39e760c5ca3b1cfa087927580-0'
        : 'biocontainers/mulled-v2-c80ae8d0fe5685926c9bc673e400ff09a71844fd:e01414b01cd5729e641b84adaf1ce4fd6181bcb8-2'}"

    input:
    tuple val(meta), path(bam)
    path bed

    output:
    tuple val(meta), path("*distance.txt"), emit: distance
    tuple val(meta), path("*freq.txt"), emit: freq
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val(meta), path("*.r"), emit: rscript
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    inner_distance.py \\
        -i ${bam} \\
        ${args} \\
        -o ${prefix} \\
        -r ${bed}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$( inner_distance.py --version | sed 's/.*.py //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.inner_distance.txt
    touch ${prefix}.inner_distance_freq.txt
    touch ${prefix}.inner_distance_plot.r
    touch ${prefix}.inner_distance_plot.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$( inner_distance.py --version | sed 's/.*.py //' )
    END_VERSIONS
    """
}
