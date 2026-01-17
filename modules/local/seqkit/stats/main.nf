#!/usr/bin/env nextflow

process SEQKIT_STATS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'
        : 'biocontainers/seqkit:2.12.0--he881be0_1'}"

    input:
    tuple val(meta), path(reads)

    output:
    env ("AVG_READLENGTH"), emit: avg_length
    path ("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def fastq_in = meta.single_end ? "${reads}" : "${reads[0]}"
    """
    AVG_READLENGTH=\$(seqkit stats ${fastq_in} | awk 'NR==2 {printf "%d\\n", \$8+0.5}')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """

    stub:
    """
    AVG_READLENGTH=100

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gxf2bed: \$( gxf2bed --version | grep 'gxf2bed [0-9]' | sed 's/.*ed //' )
    END_VERSIONS
    """
}
