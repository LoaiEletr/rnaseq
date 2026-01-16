#!/usr/bin/env nextflow

process KALLISTO_INDEX {
    tag "${transcriptome.baseName}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/kallisto:0.51.1--h2b92561_2'
        : 'biocontainers/kallisto:0.51.1--h2b92561_2'}"

    input:
    path transcriptome

    output:
    path "*.idx", emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.argg ?: ''
    def prefix = task.ext.prefix ?: "${transcriptome.baseName}"
    """
    kallisto \\
        index \\
        -t ${task.cpus} \\
        ${args} \\
        -i ${prefix}.idx \\
        ${transcriptome}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$( kallisto version | sed 's/.*ion //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${transcriptome.baseName}"
    """
    touch ${prefix}.idx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$( kallisto version | sed 's/.*ion //' )
    END_VERSIONS
    """
}
