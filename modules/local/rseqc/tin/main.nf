#!/usr/bin/env nextflow

process RSEQC_TIN {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-c80ae8d0fe5685926c9bc673e400ff09a71844fd:29c8e89bc12d33b39e760c5ca3b1cfa087927580-0'
        : 'biocontainers/mulled-v2-c80ae8d0fe5685926c9bc673e400ff09a71844fd:e01414b01cd5729e641b84adaf1ce4fd6181bcb8-2'}"

    input:
    tuple val(meta), path(bam), path(bai)
    path bed

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val(meta), path("*.xls"), emit: xls
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    tin.py \\
        -r ${bed} \\
        ${args} \\
        -i ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$( tin.py --version | sed 's/.*.py //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${bam.baseName}.tin.xls
    touch ${bam.baseName}.summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$( read_duplication.py --version | sed 's/.*.py //' )
    END_VERSIONS
    """
}
