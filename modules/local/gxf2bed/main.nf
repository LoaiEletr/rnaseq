#!/usr/bin/env nextflow

process GXF2BED {
    tag "${gtf.baseName}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gxf2bed:0.2.4--ha6fb395_0'
        : 'biocontainers/gxf2bed:0.2.7--ha6fb395_0'}"

    input:
    path gtf

    output:
    path "*.bed", emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${gtf.baseName}"
    """
    gxf2bed \\
        ${args} \\
        -i ${gtf} \\
        -o ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gxf2bed: \$( gxf2bed --version | grep 'gxf2bed [0-9]' | sed 's/.*ed //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${gtf.baseName}"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gxf2bed: \$( gxf2bed --version | grep 'gxf2bed [0-9]' | sed 's/.*ed //' )
    END_VERSIONS
    """
}
