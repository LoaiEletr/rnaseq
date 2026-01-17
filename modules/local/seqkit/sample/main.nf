#!/usr/bin/env nextflow

process SEQKIT_SAMPLE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'
        : 'biocontainers/seqkit:2.12.0--he881be0_1'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_1M.fastq.gz"), emit: reads
    tuple val(meta), path("*.log"), emit: log
    path ("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        seqkit \\
            sample \\
            -s100 \\
            ${reads} \\
            -n 1000000 \\
            ${args} \\
            -o ${prefix}_1M.fastq.gz \\
            2>| >(tee ${prefix}_1_1M.log >&2 )

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$( seqkit version | sed 's/seqkit v//' )
        END_VERSIONS
        """
    }
    else {
        """
        seqkit \\
            sample \\
            -s100 \\
            ${reads[0]} \\
            -n 1000000 \\
            ${args} \\
            -o ${prefix}_1_1M.fastq.gz \\
            2>| >(tee ${prefix}_1_1M.log >&2 )

        seqkit \\
            sample \\
            -s100 \\
            ${reads[1]} \\
            -n 1000000 \\
            -j ${task.cpus} \\
            ${args} \\
            -o ${prefix}_2_1M.fastq.gz \\
            2>| >(tee ${prefix}_2_1M.log >&2 )

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$( seqkit version | sed 's/seqkit v//' )
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        echo "" | gzip > ${prefix}_1M.fastq.gz
        touch ${prefix}_1M.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$( seqkit version | sed 's/seqkit v//' )
        END_VERSIONS
        """
    }
    else {
        """
        echo "" | gzip > ${prefix}_1_1M.fastq.gz
        echo "" | gzip > ${prefix}_2_1M.fastq.gz
        touch ${prefix}_{1,2}_1M.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$( seqkit version | sed 's/seqkit v//' )
        END_VERSIONS
        """
    }
}
