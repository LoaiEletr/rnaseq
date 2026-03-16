process SORTMERNA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2b/2bb6a75d1ac1f321bae03ad8ec6b6ce9e2303516c016d2c7d758da1b5a65f963/data'
        : 'community.wave.seqera.io/library/sortmerna:4.3.6--4946a3f74bc12d4d'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(reference_dir)
    val rrna_db_type
    val only_build_index

    output:
    tuple val(meta), path("*non_rRNA*.fq.gz"), optional: true, emit: reads
    tuple val(meta), path("index"), optional: true, emit: index
    tuple val(meta), path("*.log"), optional: true, emit: log
    tuple val("${task.process}"), val('sortmerna'), eval("sortmerna --version | grep Sort | sed 's/.*ion //'"), topic: versions, emit: versions_sortmerna

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: (only_build_index == true ? "index" : "${meta.id}")
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
            --ref ${reference_dir}/${ref_file} \\
            ${args} \\
            --workdir ./ \\
            --index 1 \\
            --idx-dir ${prefix} \\
            -v
        """
    }
    else {
        """
        sortmerna \\
            --ref ${reference_dir}/${ref_file} \\
            ${args} \\
            --workdir ./ \\
            ${fastq_in} \\
            --idx-dir ${index} \\
            --index 0 \\
            --fastx \\
            --aligned ${prefix}_rRNA_reads \\
            --other ${prefix}_non_rRNA_reads
        """
    }

    stub:
    def prefix = task.ext.prefix ?: (only_build_index == true ? "index" : "${meta.id}")
    def mkdir_index = only_build_index == true ? "mkdir -p index" : ''
    def fastq_log_out = only_build_index == false ? (meta.single_end ? "echo '' | gzip > ${prefix}_non_rRNA_reads.fq.gz ; touch rRNA_reads.log" : "echo '' | gzip > ${prefix}_non_rRNA_reads_fwd.fq.gz ; echo '' | gzip > ${prefix}_non_rRNA_reads_rev.fq.gz ; touch rRNA_reads.log") : ''
    """
    ${mkdir_index}
    ${fastq_log_out}
    """
}
