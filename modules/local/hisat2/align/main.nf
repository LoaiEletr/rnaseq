process HISAT2_ALIGN {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/16/169137662076ae5901d55f7655100f329af2e5301a0b6e77d3f80c079f155335/data'
        : 'community.wave.seqera.io/library/hisat2_samtools:3c85db656ffe50a2'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(fasta)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.summary"), emit: summary
    tuple val("${task.process}"), val('hisat2'), eval("hisat2 --version | head -n 1 | awk '{print \$3}' | sed 's/^.*version //'"), topic: versions, emit: versions_hisat2
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | head -n 1 | sed 's/samtools //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fastq_in = meta.single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    strand_flag = ""
    if (meta.lib_type in ["IU", "U", "MU", "OU"]) {
        strand_flag = ""
    }
    else if (meta.lib_type in ["SR", "ISR", "MSR", "OSR", "reverse"]) {
        strand_flag = "--rna-strandness RF"
    }
    else if (meta.lib_type in ["SF", "ISF", "MSF", "OSF", "forward"]) {
        strand_flag = "--rna-strandness FR"
    }
    else {
        error("Unknown strandness type: ${meta.lib_type}")
    }
    """
    hisat2 \\
        -x ${index}/${meta3.id} \\
        ${fastq_in} \\
        ${strand_flag} \\
        ${args} \\
        --rg-id ${meta.id} \\
        --rg "SM:${meta.id}" \\
        --rg "PL:Illumina" \\
        --new-summary \\
        --summary-file ${prefix}.hisat2.summary \\
        | samtools view -bS - > ${prefix}.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.hisat2.summary
    touch ${prefix}.bam
    """
}
