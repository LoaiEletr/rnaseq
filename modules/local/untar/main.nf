process UNTAR {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(archive)

    output:
    tuple val(meta), path("${prefix}"), emit: untar
    tuple val("${task.process}"), val('untar'), eval("tar --version | head -n 1 | sed 's/.*) //'"), topic: versions, emit: versions_untar

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = archive.toString() - '.tar.gz'
    """
    tar \\
        -xzvf \\
        ${args} \\
        ${archive} \\
        --one-top-level=${prefix}
    """

    stub:
    prefix = archive.toString() - '.tar.gz'
    """
    touch ${prefix}
    """
}
