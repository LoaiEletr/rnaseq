process MULTIQC {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/34/34e733a9ae16a27e80fe00f863ea1479c96416017f24a907996126283e7ecd4d/data'
        : 'community.wave.seqera.io/library/multiqc:1.33--ee7739d47738383b'}"

    input:
    tuple val(meta), path(multiqc_files)
    tuple val(meta2), path(multiqc_config)

    output:
    tuple val(meta), path("*.html"), emit: report
    tuple val(meta), path("*_data"), emit: data
    tuple val(meta), path("*_plots"), optional: true, emit: plots
    tuple val("${task.process}"), val('multiqc'), eval("multiqc --version | sed 's/.*version //'"), emit: versions_multiqc

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def config = multiqc_config ? "--config ${multiqc_config}" : ''
    """
    multiqc \\
        --force \\
        ${args} \\
        ${config} \\
        ${multiqc_files} \\
    """

    stub:
    """
    mkdir multiqc_data
    mkdir multiqc_plots
    touch multiqc_report.html
    """
}
