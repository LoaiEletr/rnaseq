#!/usr/bin/env nextflow

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
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_string_ppi.R ${wgcna_results_rds} ${species_name} ${score_threshold}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-stringdb: \$(Rscript -e "library(STRINGdb); cat(as.character(packageVersion('STRINGdb')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir STRING_PPI

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-stringdb: \$(Rscript -e "library(STRINGdb); cat(as.character(packageVersion('STRINGdb')))")
    END_VERSIONS
    """
}
