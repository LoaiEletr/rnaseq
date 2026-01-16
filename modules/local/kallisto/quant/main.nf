#!/usr/bin/env nextflow

process KALLISTO_QUANT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/kallisto:0.51.1--h2b92561_2'
        : 'biocontainers/kallisto:0.51.1--h2b92561_2'}"

    input:
    tuple val(meta), path(reads)
    path index
    path gtf
    val bootstrap_count
    val fragment_length
    val fragment_length_sd

    output:
    tuple val(meta), path("${prefix}"), emit: quant_dir
    tuple val(meta), path("${prefix}/run_info.json"), emit: json_info
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def bootstrap_count_input = bootstrap_count ? "-b ${bootstrap_count}" : ''
    def gtf_input = gtf ? "--gtf ${gtf}" : ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def fastq_in = meta.single_end ? "--single -l ${fragment_length} -s ${fragment_length_sd} ${reads}" : "${reads[0]} ${reads[1]}"
    strand_flag = ""

    // Use the meta.lib_type input value (from auto-detection)
    if (meta.lib_type == "IU" || meta.lib_type == "U") {
        strand_flag = ""
    }
    else if (meta.lib_type in ["ISR", "SF", "MSR", "OSR", "reverse"]) {
        strand_flag = "--rf-stranded"
    }
    else if (meta.lib_type in ["ISF", "SF", "MSF", "OSF", "forward"]) {
        strand_flag = "--fr-stranded"
    }
    else {
        error("Unknown meta.lib_type type: ${meta.lib_type}")
    }
    """
    kallisto \\
        quant \\
        ${args} \\
        ${bootstrap_count_input} \\
        ${gtf_input} \\
        ${strand_flag} \\
        -t ${task.cpus} \\
        -i ${index} \\
        -o ${prefix} \\
        ${fastq_in} \\
        2>| >( tee kallisto_${prefix}_quant.log >&2 )

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$( kallisto version | sed 's/.*ion //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch kallisto_${prefix}_quant.log
    touch ${prefix}/abundance.h5
    touch ${prefix}/abundance.tsv
    touch ${prefix}/run_info.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$( kallisto version | sed 's/.*ion //' )
    END_VERSIONS
    """
}
