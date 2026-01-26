#!/usr/bin/env nextflow

process ISOFORMSWITCHANALYZER {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3e/3e73e76384587e47c557b88bfce8212a821888d52b4647967e69895fe64f2083/data'
        : 'community.wave.seqera.io/library/bioconductor-isoformswitchanalyzer_bioconductor-rhdf5_r-base_r-ggplot2_r-mass:e41cc44dc5deb629'}"

    input:
    path quant_dir
    path samplesheet
    val quant_type
    path gtf
    path transcript
    val method
    val dif_cutoff
    val pval_threshold
    val ntop_isoforms

    output:
    path "isoform_output", emit: isoform_results, optional: true
    path "isoform_output/significant_DIU_gene_ids.rds", emit: significant_geneids, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_isoform_switch_analysis.R ${quant_dir} ${samplesheet} ${quant_type} ${gtf} ${transcript} ${method} ${dif_cutoff} ${pval_threshold} ${ntop_isoforms}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-isoformswitchanalyzer: \$(Rscript -e "library(IsoformSwitchAnalyzeR); cat(as.character(packageVersion('IsoformSwitchAnalyzeR')))")
        bioconductor-rhdf5: \$(Rscript -e "library(rhdf5); cat(as.character(packageVersion('rhdf5')))")
        r-mass: \$(Rscript -e "library(MASS); cat(as.character(packageVersion('MASS')))")
        r-ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir isoform_output
    touch isoform_output/significant_DIU_gene_ids.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-isoformswitchanalyzer: \$(Rscript -e "library(IsoformSwitchAnalyzeR); cat(as.character(packageVersion('IsoformSwitchAnalyzeR')))")
        bioconductor-rhdf5: \$(Rscript -e "library(rhdf5); cat(as.character(packageVersion('rhdf5')))")
        r-mass: \$(Rscript -e "library(MASS); cat(as.character(packageVersion('MASS')))")
        r-ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
    END_VERSIONS
    """
}
