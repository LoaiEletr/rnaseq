process STRING_PPI {
    tag "${wgcna_results_rds}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9f/9f256b97adf0f7548a53d99157868eccc6730478068819f13905327a6fac1c87/data'
        : 'community.wave.seqera.io/library/bioconductor-stringdb_r-base:6296d2774720661c'}"

    input:
    path wgcna_results_rds
    val species_name
    val score_threshold

    output:
    path "STRING_PPI", emit: string_results
    tuple val("${task.process}"), val('r-base'), eval("Rscript -e 'cat(as.character(getRversion()))'"), topic: versions, emit: versions_rbase
    tuple val("${task.process}"), val('bioconductor-stringdb'), eval("Rscript -e \"cat(as.character(packageVersion('STRINGdb')))\""), topic: versions, emit: versions_stringdb

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_string_ppi.R ${wgcna_results_rds} ${species_name} ${score_threshold}
    """

    stub:
    """
    mkdir STRING_PPI
    """
}
