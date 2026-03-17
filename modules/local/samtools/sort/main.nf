process SAMTOOLS_SORT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b9/b946e2f0e77ec69853787dfc8b312bd7e9d5c65a11a613ce918469a9566992e3/data'
        : 'community.wave.seqera.io/library/samtools:1.23--12d9384dd0649f36'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.sorted.bam"), emit: sorted_bam
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | head -n 1 | sed 's/samtools //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: (bam.name.contains("umi_dedup") ? "${meta.id}.umi_dedup" : "${meta.id}")
    """
    samtools \\
        sort ${bam} \\
        -o ${prefix}.sorted.bam \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: (bam.name.contains("umi_dedup") ? "${meta.id}.umi_dedup" : "${meta.id}")
    """
    touch ${prefix}.sorted.bam
    """
}
