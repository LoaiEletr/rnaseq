process GUNZIP {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(archive)

    output:
    tuple val(meta), path("${gunzip}"), emit: gunzip
    tuple val("${task.process}"), val('gunzip'), eval("gunzip --version | grep 'gunzip' | awk '{print \$3}' | sed 's/^.*(gzip) //'"), topic: versions, emit: versions_gunzip

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    gunzip = archive.toString() - '.gz'
    """
    gzip \\
        -cd \\
        ${args} \\
        ${archive} \\
        > ${gunzip}
    """

    stub:
    gunzip = archive.toString() - '.gz'
    """
    touch ${gunzip}
    """
}
