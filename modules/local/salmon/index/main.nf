#!/usr/bin/env nextflow

process SALMON_INDEX {
    tag "${transcriptome}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/salmon:1.10.3--haf24da9_3'
        : 'biocontainers/salmon:1.10.3--h45fbf2d_5'}"

    input:
    path transcriptome

    output:
    path ("transcripts_index"), emit: index
    path ("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "transcripts_index"
    """
    salmon \\
        index \\
        -t ${transcriptome} \\
        -p ${task.cpus} \\
        -i ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$( salmon --version | sed 's/salmon //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "transcripts_index"
    """
    mkdir -p ${prefix}
    touch ${prefix}/complete_ref_lens.bin
    touch ${prefix}/ctable.bin
    touch ${prefix}/ctg_offsets.bin
    touch ${prefix}/duplicate_clusters.tsv
    touch ${prefix}/info.json
    touch ${prefix}/mphf.bin
    touch ${prefix}/pos.bin
    touch ${prefix}/pre_indexing.log
    touch ${prefix}/rank.bin
    touch ${prefix}/refAccumLengths.bin
    touch ${prefix}/ref_indexing.log
    touch ${prefix}/reflengths.bin
    touch ${prefix}/refseq.bin
    touch ${prefix}/seq.bin
    touch ${prefix}/versionInfo.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$( salmon --version | sed 's/salmon //' )
    END_VERSIONS
    """
}
