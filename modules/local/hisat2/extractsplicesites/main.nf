process HISAT2_EXTRACTSPLICESITES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/eb/eb454805f1aad64254ef3d4a806ce8ffc75605726715ca7328e2db224903b700/data'
        : 'community.wave.seqera.io/library/hisat2:2.2.1--df34d2bb25ac6de5'}"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.splice_sites.tsv"), emit: splice_sites
    tuple val("${task.process}"), val('hisat2'), eval("hisat2 --version | head -n 1 | awk '{print \$3}' | sed 's/^.*version //'"), topic: versions, emit: versions_hisat2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hisat2_extract_splice_sites.py \\
        ${gtf} \\
        ${args} \\
        > ${prefix}.splice_sites.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.splice_sites.tsv
    """
}
