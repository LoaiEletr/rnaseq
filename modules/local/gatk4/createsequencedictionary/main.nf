#!/usr/bin/env nextflow

process GATK4_CREATESEQUENCEDICTIONARY {
    tag "${fasta}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    path fasta

    output:
    path "*.dict", emit: dict
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def avail_mem = task.memory.mega
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        CreateSequenceDictionary \\
        ${args} \\
        --REFERENCE ${fasta} \\
        --URI ${fasta} \\
        --TMP_DIR .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$( gatk --version | sed -n '/GATK.*v/s/.*v//p' )
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta.baseName}.dict

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$( gatk --version | sed -n '/GATK.*v/s/.*v//p' )
    END_VERSIONS
    """
}
