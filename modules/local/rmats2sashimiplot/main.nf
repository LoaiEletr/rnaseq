#!/usr/bin/env nextflow

process RMATS2SASHIMIPLOT {
    tag "${event_type}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/rmats2sashimiplot:3.0.0--py39hdff8610_2'
        : 'biocontainers/rmats2sashimiplot:3.0.0--py310ha6fa2df_2'}"

    input:
    tuple val(event_type), path(rmats_sigevent), path(bam_list1), path(bam_list2)
    path bam_files

    output:
    path ("${prefix}"), emit: sashimi_output
    path ("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "sashimi_${event_type}_out"
    if (!(event_type in ["SE", "MXE", "A3SS", "A5SS", "RI"])) {
        error("Unknown event type: ${prefix}")
    }
    def condition_flag = bam_list1.baseName.toLowerCase() =~ /control/ ? "--l1 ${bam_list2.baseName} --l2 ${bam_list1.baseName}" : "--l1 ${bam_list1.baseName} --l2 ${bam_list2.baseName}"
    def bam_list_flag = bam_list1.baseName.toLowerCase() =~ /control/ ? "--b1 ${bam_list2} --b2 ${bam_list1}" : "--b1 ${bam_list1} --b2 ${bam_list2}"
    """
    rmats2sashimiplot \\
        -o ${prefix} \\
        ${condition_flag} \\
        --event-type ${event_type} \\
        -e ${rmats_sigevent} \\
        ${bam_list_flag}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rmats2sashimiplot: \$( pip show rmats2sashimiplot | grep "Version:" | awk '{print \$2}' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "sashimi_${event_type}_out"
    """
    mkdir ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rmats2sashimiplot: \$( pip show rmats2sashimiplot | grep "Version:" | awk '{print \$2}' )
    END_VERSIONS
    """
}
