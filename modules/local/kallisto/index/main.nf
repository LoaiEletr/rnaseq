process KALLISTO_INDEX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/70/70a7a1b2f5dc1ca0e801e958ab14ea029f662f1ead026ba9cdff59f99995f19c/data'
        : 'community.wave.seqera.io/library/kallisto:0.51.1--b63691b6841c7a52'}"

    input:
    tuple val(meta), path(transcriptome)

    output:
    tuple val(meta), path("*.idx"), emit: index
    tuple val("${task.process}"), val('kallisto'), eval("kallisto version | sed 's/.*ion //'"), topic: versions, emit: versions_kallisto

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.argg ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    kallisto \\
        index \\
        -t ${task.cpus} \\
        ${args} \\
        -i ${prefix}.idx \\
        ${transcriptome}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.idx
    """
}
