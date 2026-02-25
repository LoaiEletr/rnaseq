process MAKE_BAMLIST {
    tag "${meta.condition}"
    label 'process_low'

    input:
    tuple val(meta), path(bam_files)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val("${task.process}"), val('bash'), eval("bash --version | head -n1 | awk '{print \$4}' | sed 's/^.*version //; s/(.*//'"), topic: versions, emit: versions_bash

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.condition}"
    """
    echo \\
        ${bam_files.join(',')} \\
        > ${prefix}.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.condition}"
    """
    touch ${prefix}.txt
    """
}
