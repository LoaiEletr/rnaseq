#!/usr/bin/env nextflow

process KEGG_ENRICHMENT {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b2b1b22d2d1fd80ad2a9ce1feb69510340ab48504e1de13832bbbf3c12166fe0/data'
        : 'community.wave.seqera.io/library/bioconductor-biomart_bioconductor-clusterprofiler_bioconductor-org.at.tair.db_bioconductor-org.bt.eg.db_pruned:263fa5783ae58d65'}"

    input:
    path unfiltered_deg
    val pvalue_threshold
    val logfc_threshold
    val species_name
    path wgcna_list
    path masigpro_list
    path genes_diu
    path genes_deu
    val ntop_processes

    output:
    path "KEGG_results", optional: true, emit: kegg_results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def unfiltered_deg_flag = unfiltered_deg ?: null
    def wgcna_list_flag = wgcna_list ?: null
    def masigpro_list_flag = masigpro_list ?: null
    def logfc_threshold_flag = logfc_threshold ?: null
    def genes_diu_flag = genes_diu ?: null
    def genes_deu_flag = genes_deu ?: null
    """
    run_kegg_enrichment.R ${unfiltered_deg_flag} ${pvalue_threshold} ${logfc_threshold_flag} ${species_name} ${wgcna_list_flag} ${masigpro_list_flag} ${genes_diu_flag} ${genes_deu_flag} ${ntop_processes}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-clusterprofiler: \$(Rscript -e "library(clusterProfiler); cat(as.character(packageVersion('clusterProfiler')))")
        r-ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        bioconductor-org.hs.eg.db: \$(Rscript -e "library(org.Hs.eg.db); cat(as.character(packageVersion('org.Hs.eg.db')))")
        r-purrr: \$(Rscript -e "library(purrr); cat(as.character(packageVersion('purrr')))")
        bioconductor-biomart: \$(Rscript -e "library(biomaRt); cat(as.character(packageVersion('biomaRt')))")
        bioconductor-org.mm.eg.db: \$(Rscript -e "library(org.Mm.eg.db); cat(as.character(packageVersion('org.Mm.eg.db')))")
        bioconductor-org.rn.eg.db: \$(Rscript -e "library(org.Rn.eg.db); cat(as.character(packageVersion('org.Rn.eg.db')))")
        r-dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        r-tibble: \$(Rscript -e "library(tibble); cat(as.character(packageVersion('tibble')))")
        bioconductor-org.sc.sgd.db: \$(Rscript -e "library(org.Sc.sgd.db); cat(as.character(packageVersion('org.Sc.sgd.db')))")
        bioconductor-org.dm.eg.db: \$(Rscript -e "library(org.Dm.eg.db); cat(as.character(packageVersion('org.Dm.eg.db')))")
        bioconductor-org.dr.eg.db: \$(Rscript -e "library(org.Dr.eg.db); cat(as.character(packageVersion('org.Dr.eg.db')))")
        bioconductor-org.ce.eg.db: \$(Rscript -e "library(org.Ce.eg.db); cat(as.character(packageVersion('org.Ce.eg.db')))")
        bioconductor-org.at.tair.db: \$(Rscript -e "library(org.At.tair.db); cat(as.character(packageVersion('org.At.tair.db')))")
        bioconductor-org.gg.eg.db: \$(Rscript -e "library(org.Gg.eg.db); cat(as.character(packageVersion('org.Gg.eg.db')))")
        bioconductor-org.bt.eg.db: \$(Rscript -e "library(org.Bt.eg.db); cat(as.character(packageVersion('org.Bt.eg.db')))")
        bioconductor-org.ss.eg.db: \$(Rscript -e "library(org.Ss.eg.db); cat(as.character(packageVersion('org.Ss.eg.db')))")
        bioconductor-org.cf.eg.db: \$(Rscript -e "library(org.Cf.eg.db); cat(as.character(packageVersion('org.Cf.eg.db')))")
        bioconductor-org.mmu.eg.db: \$(Rscript -e "library(org.Mmu.eg.db); cat(as.character(packageVersion('org.Mmu.eg.db')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir KEGG_results

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-clusterprofiler: \$(Rscript -e "library(clusterProfiler); cat(as.character(packageVersion('clusterProfiler')))")
        r-ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        bioconductor-org.hs.eg.db: \$(Rscript -e "library(org.Hs.eg.db); cat(as.character(packageVersion('org.Hs.eg.db')))")
        r-purrr: \$(Rscript -e "library(purrr); cat(as.character(packageVersion('purrr')))")
        bioconductor-biomart: \$(Rscript -e "library(biomaRt); cat(as.character(packageVersion('biomaRt')))")
        bioconductor-org.mm.eg.db: \$(Rscript -e "library(org.Mm.eg.db); cat(as.character(packageVersion('org.Mm.eg.db')))")
        bioconductor-org.rn.eg.db: \$(Rscript -e "library(org.Rn.eg.db); cat(as.character(packageVersion('org.Rn.eg.db')))")
        r-dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        r-tibble: \$(Rscript -e "library(tibble); cat(as.character(packageVersion('tibble')))")
        bioconductor-org.sc.sgd.db: \$(Rscript -e "library(org.Sc.sgd.db); cat(as.character(packageVersion('org.Sc.sgd.db')))")
        bioconductor-org.dm.eg.db: \$(Rscript -e "library(org.Dm.eg.db); cat(as.character(packageVersion('org.Dm.eg.db')))")
        bioconductor-org.dr.eg.db: \$(Rscript -e "library(org.Dr.eg.db); cat(as.character(packageVersion('org.Dr.eg.db')))")
        bioconductor-org.ce.eg.db: \$(Rscript -e "library(org.Ce.eg.db); cat(as.character(packageVersion('org.Ce.eg.db')))")
        bioconductor-org.at.tair.db: \$(Rscript -e "library(org.At.tair.db); cat(as.character(packageVersion('org.At.tair.db')))")
        bioconductor-org.gg.eg.db: \$(Rscript -e "library(org.Gg.eg.db); cat(as.character(packageVersion('org.Gg.eg.db')))")
        bioconductor-org.bt.eg.db: \$(Rscript -e "library(org.Bt.eg.db); cat(as.character(packageVersion('org.Bt.eg.db')))")
        bioconductor-org.ss.eg.db: \$(Rscript -e "library(org.Ss.eg.db); cat(as.character(packageVersion('org.Ss.eg.db')))")
        bioconductor-org.cf.eg.db: \$(Rscript -e "library(org.Cf.eg.db); cat(as.character(packageVersion('org.Cf.eg.db')))")
        bioconductor-org.mmu.eg.db: \$(Rscript -e "library(org.Mmu.eg.db); cat(as.character(packageVersion('org.Mmu.eg.db')))")
    END_VERSIONS
    """
}
