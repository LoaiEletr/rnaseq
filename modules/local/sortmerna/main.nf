#!/usr/bin/env nextflow

process SORTMERNA {
    tag "${prefix}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sortmerna:4.3.4--h9ee0642_0'
        : 'community.wave.seqera.io/library/sortmerna:4.3.6--4946a3f74bc12d4d'}"

    input:
    tuple val(meta), path(reads)
    path index
    path ref_dir
    val rrna_db_type
    val only_build_index

    output:
    tuple val(meta), path("*non_rRNA*.fq.gz"), optional: true, emit: reads
    path ("index"), optional: true, emit: index
    tuple val(meta), path("*.log"), optional: true, emit: log
    path ("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: (only_build_index == true ? "index" : "${meta.id}")
    def args = task.ext.args ?: ''
    def fastq_in = meta.single_end ? "--reads ${reads}" : "--paired_in --out2 --reads ${reads[0]} --reads ${reads[1]}"
    def ref_file = (rrna_db_type == "fast"
        ? "smr_v4.3_fast_db.fasta"
        : rrna_db_type == "default"
            ? "smr_v4.3_default_db.fasta"
            : rrna_db_type == "sensitive"
                ? "smr_v4.3_sensitive_db.fasta"
                : "smr_v4.3_sensitive_db_rfam_seeds.fasta")
    if (only_build_index) {
        """
        sortmerna \\
            --ref ${ref_dir}/${ref_file} \\
            ${args} \\
            --index 1 \\
            --idx-dir ${prefix} \\
            -v

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sortmerna: \$( sortmerna --version | grep Sort | sed 's/.*ion //' )
        END_VERSIONS
        """
    }
    else {
        """
        sortmerna \\
            --ref ${ref_dir}/${ref_file} \\
            ${args} \\
            --workdir ./ \\
            ${fastq_in} \\
            --idx-dir ${index} \\
            --index 0 \\
            --fastx \\
            --aligned ${prefix}_rRNA_reads \\
            --other ${prefix}_non_rRNA_reads

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sortmerna: \$( sortmerna --version | grep Sort | sed 's/.*ion //' )
        END_VERSIONS
        """
    }

    stub:
    prefix = task.ext.prefix ?: (only_build_index == true ? "index" : "${meta.id}")
    def mkdir_index = only_build_index == true ? "mkdir -p index" : ''
    def fastq_log_out = only_build_index == false ? (meta.single_end ? "echo '' | gzip > ${prefix}_non_rRNA_reads.fq.gz ; touch rRNA_reads.log" : "echo '' | gzip > ${prefix}_non_rRNA_reads_fwd.fq.gz ; echo '' | gzip > ${prefix}_non_rRNA_reads_rev.fq.gz ; touch rRNA_reads.log") : ''
    """
    ${mkdir_index}
    ${fastq_log_out}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sortmerna: \$( sortmerna --version | grep Sort | sed 's/.*ion //' )
    END_VERSIONS
    """
}
