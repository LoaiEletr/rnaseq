#!/usr/bin/env nextflow

process HISAT2_ALIGN {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:0e773bb207600fcb4d38202226eb20a33c7909b6-0'
        : 'biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0'}"

    input:
    tuple val(meta), path(reads)
    path index
    path ref

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.summary"), emit: summary
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fastq_in = meta.single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    strand_flag = ""
    if (meta.lib_type in ["IU", "U", "MU", "OU"]) {
        strand_flag = ""
    }
    else if (meta.lib_type in ["SR", "ISR", "MSR", "OSR", "reverse"]) {
        strand_flag = "--rna-strandness RF"
    }
    else if (meta.lib_type in ["SF", "ISF", "MSF", "OSF", "forward"]) {
        strand_flag = "--rna-strandness FR"
    }
    else {
        error("Unknown strandness type: ${meta.lib_type}")
    }
    """
    hisat2 \\
        -x ${index}/${ref.getSimpleName()} \\
        ${fastq_in} \\
        ${strand_flag} \\
        ${args} \\
        --summary-file ${prefix}.hisat2.summary \\
        | samtools view -bS - > ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$( hisat2 --version | head -n 1 | awk '{print \$3}' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.hisat2.summary
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$( hisat2 --version | head -n 1 | awk '{print \$3}' )
    END_VERSIONS
    """
}
