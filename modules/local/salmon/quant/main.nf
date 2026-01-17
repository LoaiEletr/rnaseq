#!/usr/bin/env nextflow

process SALMON_QUANT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/salmon:1.10.3--haf24da9_3'
        : 'biocontainers/salmon:1.10.3--h45fbf2d_5'}"

    input:
    tuple val(meta), path(reads)
    path index
    path gtf
    val libtype

    output:
    tuple val(meta), path("${prefix}"), emit: quant_dir, optional: true
    tuple val(meta), env("SALMON_STRANDNESS"), emit: strandness, optional: true
    tuple val(meta), path("*.log"), emit: log
    path ("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def fastq_in = meta.single_end ? "-r ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"

    def paired_end_strandedness_opts = [
        'A',
        'IS',
        'IU',
        'ISF',
        'ISR',
        'OS',
        'OU',
        'OSF',
        'OSR',
        'MS',
        'MU',
        'MSF',
        'MSR',
    ]

    def single_end_strandedness_opts = [
        'U',
        'A',
        'ISF',
        'ISR',
    ]

    // Determine strand flag
    if (meta.single_end) {
        // Single-end library
        if (libtype) {
            if (single_end_strandedness_opts.contains(libtype)) {
                strand_flag = "-l ${libtype}"
            }
            else {
                error("❌ ERROR: Unknown library type '${libtype}' provided for single-end reads. Allowed: ${single_end_strandedness_opts.join(', ')}")
            }
        }
        else {
            if (meta.lib_type.toLowerCase() == "forward") {
                strand_flag = "-l SF"
            }
            else if (meta.lib_type.toLowerCase() == "reverse") {
                strand_flag = "-l SR"
            }
            else if (single_end_strandedness_opts.contains(meta.lib_type)) {
                strand_flag = "-l ${meta.lib_type}"
            }
            else {
                error("Unknown meta.lib_type provided for single-end: ${meta.lib_type}. Use one of: ${single_end_strandedness_opts.join(', ')}")
            }
        }
    }
    else {
        // Paired-end library
        if (libtype) {
            if (paired_end_strandedness_opts.contains(libtype)) {
                strand_flag = "-l ${libtype}"
            }
            else {
                error("❌ ERROR: Unknown library type '${libtype}' provided for paired-end reads. Allowed: ${paired_end_strandedness_opts.join(', ')}")
            }
        }
        else {
            if (meta.lib_type.toLowerCase() == "forward") {
                strand_flag = "-l ISF"
            }
            else if (meta.lib_type.toLowerCase() == "reverse") {
                strand_flag = "-l ISR"
            }
            else if (paired_end_strandedness_opts.contains(meta.lib_type)) {
                strand_flag = "-l ${meta.lib_type}"
            }
            else {
                error("Unknown meta.lib_type provided for paired-end: ${meta.lib_type}. Use one of: ${paired_end_strandedness_opts.join(', ')}")
            }
        }
    }

    """
    salmon quant \\
        -i ${index} \\
        ${args} \\
        ${strand_flag} \\
        --geneMap ${gtf} \\
        -p ${task.cpus} \\
        ${fastq_in} \\
        -o ${prefix} \\
        2>| >( tee ${prefix}_quant.log >&2 )

    SALMON_STRANDNESS=\$(grep "most likely library type as" ${prefix}_quant.log | awk -F'as ' '{print \$2}' || echo "Not detected")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$( salmon --version | sed 's/salmon //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix} ${prefix}/aux_info ${prefix}/libParams ${prefix}/logs

    touch ${prefix}_quant.log
    touch ${prefix}/cmd_info.json
    touch ${prefix}/lib_format_counts.json
    touch ${prefix}/quant.sf

    touch ${prefix}/aux_info/ambig_info.tsv
    touch ${prefix}/aux_info/meta_info.json

    echo "" | gzip > ${prefix}/aux_info/expected_bias.gz
    echo "" | gzip > ${prefix}/aux_info/fld.gz
    echo "" | gzip > ${prefix}/aux_info/observed_bias_3p.gz
    echo "" | gzip > ${prefix}/aux_info/observed_bias.gz

    touch ${prefix}/libParams/flenDist.txt
    touch ${prefix}/logs/salmon_quant.log

    SALMON_STRANDNESS="IU"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$( salmon --version | sed 's/salmon //' )
    END_VERSIONS
    """
}
