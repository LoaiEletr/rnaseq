process HISAT2_EXTRACTSPLICESITES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/hisat2:2.2.1--h503566f_8'
        : 'biocontainers/hisat2:2.2.1--h503566f_8'}"

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
