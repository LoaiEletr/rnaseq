process MERGE_COUNTS {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/r-base:4.4.1'
        : 'biocontainers/r-base:4.4.1'}"

    input:
    path counts

    output:
    path "merged_counts.rds", emit: counts_table
    tuple val("${task.process}"), val('r-base'), eval("Rscript -e 'cat(as.character(getRversion()))'"), topic: versions, emit: versions_rbase

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    merge_counts.R ${counts}
    """

    stub:
    """
    touch merged_counts.rds
    """
}
