process FASTQC {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e0/e0c976cb2eca5fee72618a581537a4f8ea42fcae24c9b201e2e0f764fd28648a/data'
        : 'community.wave.seqera.io/library/fastqc:0.12.1--af7a5314d5015c29'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"), emit: zip
    tuple val("${task.process}"), val('fastqc'), eval("fastqc --version | sed 's/.*v//'"), topic: versions, emit: versions_fastqc

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fastqc \\
        ${args} \\
        --threads ${task.cpus} \\
        ${reads} \\
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def touch_fastqc_output = meta.single_end ? "touch ${prefix}_fastqc.html ; touch ${prefix}_fastqc.zip" : "touch ${prefix}_1_fastqc.html  ; touch ${prefix}_2_fastqc.html ; touch ${prefix}_1_fastqc.zip  ; touch ${prefix}_2_fastqc.zip"
    """
    ${touch_fastqc_output}
    """
}
