#!/usr/bin/env nextflow

process BBMAP_BBSPLIT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bbmap:39.13--he5f24ec_1'
        : 'biocontainers/bbmap:39.50--he5f24ec_0'}"

    input:
    tuple val(meta), path(reads)
    path index
    path primary_ref
    path other_ref
    val only_build_index

    output:
    tuple val(meta), path("*primary*.fastq.gz"), optional: true, emit: reads
    path "bbsplit", optional: true, emit: index
    tuple val(meta), path("*stats.txt"), optional: true, emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: (only_build_index == true ? "build_index" : "${meta.id}")
    def primary_ref_name = primary_ref ? primary_ref.getSimpleName() : ''
    def other_ref_name = other_ref ? other_ref.getSimpleName() : ''
    def args = task.ext.args ?: ''
    def avail_mem = task.memory.giga
    def fastq_in = meta.single_end ? "in=${reads}" : "in=${reads[0]} in2=${reads[1]}"
    def fastq_out = meta.single_end ? "basename=${prefix}_%.fastq.gz" : "basename=${prefix}_%_#.fastq.gz"
    if (only_build_index) {
        """
        bbsplit.sh \\
            -Xmx${avail_mem}g \\
            -Xms${avail_mem}g \\
            ${args} \\
            ref_primary_${primary_ref_name}=${primary_ref} \\
            ref_${other_ref_name}=${other_ref} \\
            path=bbsplit \\
            build=1

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bbmap: \$( bbversion.sh --version | head -n 1 )
        END_VERSIONS
        """
    }
    else {
        """
        bbsplit.sh \\
            -Xmx${avail_mem}g \\
            -Xms${avail_mem}g \\
            ${args} \\
            ${fastq_in} \\
            path=${index} \\
            ${fastq_out} \\
            maxindel=1000000 \\
            minhits=1 \\
            minratio=0.5 \\
            refstats=${prefix}.bbsplit_stats.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bbmap: \$( bbversion.sh --version | head -n 1 )
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: (only_build_index == true ? "build_index" : "${meta.id}")
    def primary_ref_name = primary_ref ? primary_ref.getSimpleName() : ''
    def gzip_fastq_out = meta.single_end ? "echo '' | gzip > ${prefix}_primary_${primary_ref_name}.fastq.gz" : "echo '' | gzip > ${prefix}_primary_${primary_ref_name}_1.fastq.gz ; echo '' | gzip > ${prefix}_primary_${primary_ref_name}_2.fastq.gz"
    if (only_build_index) {
        """
        mkdir -p bbsplit/genome/1
        mkdir -p bbsplit/index/1

        echo "" | gzip > bbsplit/genome/1/chr1.chrom.gz
        echo "" | gzip > bbsplit/genome/1/merged_ref_9219829362058850410.fa.gz
        echo "" | gzip > bbsplit/genome/1/scaffolds.txt.gz
        echo "" | gzip > bbsplit/index/1/chr1_index_k13_c5_b1.block2.gz

        touch bbsplit/genome/1/info.txt
        touch bbsplit/genome/1/namelist.txt
        touch bbsplit/genome/1/reflist.txt
        touch bbsplit/genome/1/summary.txt
        touch bbsplit/index/1/chr1_index_k13_c5_b1.block

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bbmap: \$( bbversion.sh --version | head -n 1 )
        END_VERSIONS
        """
    }
    else {
        """
        ${gzip_fastq_out}
        touch ${prefix}.bbsplit_stats.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bbmap: \$( bbversion.sh --version | head -n 1 )
        END_VERSIONS
        """
    }
}
