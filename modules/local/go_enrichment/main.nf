process GO_ENRICHMENT {
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
    path "GO_results", emit: go_results, optional: true
    tuple val("${task.process}"), val('r-base'), eval("Rscript -e 'cat(as.character(getRversion()))'"), topic: versions, emit: versions_rbase
    tuple val("${task.process}"), val('bioconductor-clusterprofiler'), eval("Rscript -e \"cat(as.character(packageVersion('clusterProfiler')))\""), topic: versions, emit: versions_clusterprofiler
    tuple val("${task.process}"), val('r-ggplot2'), eval("Rscript -e \"cat(as.character(packageVersion('ggplot2')))\""), topic: versions, emit: versions_ggplot2
    tuple val("${task.process}"), val('bioconductor-org.hs.eg.db'), eval("Rscript -e \"cat(as.character(packageVersion('org.Hs.eg.db')))\""), topic: versions, emit: versions_org_hs_eg_db
    tuple val("${task.process}"), val('r-purrr'), eval("Rscript -e \"cat(as.character(packageVersion('purrr')))\""), topic: versions, emit: versions_purrr
    tuple val("${task.process}"), val('bioconductor-biomart'), eval("Rscript -e \"cat(as.character(packageVersion('biomaRt')))\""), topic: versions, emit: versions_biomart
    tuple val("${task.process}"), val('bioconductor-org.mm.eg.db'), eval("Rscript -e \"cat(as.character(packageVersion('org.Mm.eg.db')))\""), topic: versions, emit: versions_org_mm_eg_db
    tuple val("${task.process}"), val('bioconductor-org.rn.eg.db'), eval("Rscript -e \"cat(as.character(packageVersion('org.Rn.eg.db')))\""), topic: versions, emit: versions_org_rn_eg_db
    tuple val("${task.process}"), val('r-dplyr'), eval("Rscript -e \"cat(as.character(packageVersion('dplyr')))\""), topic: versions, emit: versions_dplyr
    tuple val("${task.process}"), val('r-tibble'), eval("Rscript -e \"cat(as.character(packageVersion('tibble')))\""), topic: versions, emit: versions_tibble
    tuple val("${task.process}"), val('bioconductor-org.sc.sgd.db'), eval("Rscript -e \"cat(as.character(packageVersion('org.Sc.sgd.db')))\""), topic: versions, emit: versions_org_sc_sgd_db
    tuple val("${task.process}"), val('bioconductor-org.dm.eg.db'), eval("Rscript -e \"cat(as.character(packageVersion('org.Dm.eg.db')))\""), topic: versions, emit: versions_org_dm_eg_db
    tuple val("${task.process}"), val('bioconductor-org.dr.eg.db'), eval("Rscript -e \"cat(as.character(packageVersion('org.Dr.eg.db')))\""), topic: versions, emit: versions_org_dr_eg_db
    tuple val("${task.process}"), val('bioconductor-org.ce.eg.db'), eval("Rscript -e \"cat(as.character(packageVersion('org.Ce.eg.db')))\""), topic: versions, emit: versions_org_ce_eg_db
    tuple val("${task.process}"), val('bioconductor-org.at.tair.db'), eval("Rscript -e \"cat(as.character(packageVersion('org.At.tair.db')))\""), topic: versions, emit: versions_org_at_tair_db
    tuple val("${task.process}"), val('bioconductor-org.gg.eg.db'), eval("Rscript -e \"cat(as.character(packageVersion('org.Gg.eg.db')))\""), topic: versions, emit: versions_org_gg_eg_db
    tuple val("${task.process}"), val('bioconductor-org.bt.eg.db'), eval("Rscript -e \"cat(as.character(packageVersion('org.Bt.eg.db')))\""), topic: versions, emit: versions_org_bt_eg_db
    tuple val("${task.process}"), val('bioconductor-org.ss.eg.db'), eval("Rscript -e \"cat(as.character(packageVersion('org.Ss.eg.db')))\""), topic: versions, emit: versions_org_ss_eg_db
    tuple val("${task.process}"), val('bioconductor-org.cf.eg.db'), eval("Rscript -e \"cat(as.character(packageVersion('org.Cf.eg.db')))\""), topic: versions, emit: versions_org_cf_eg_db
    tuple val("${task.process}"), val('bioconductor-org.mmu.eg.db'), eval("Rscript -e \"cat(as.character(packageVersion('org.Mmu.eg.db')))\""), topic: versions, emit: versions_org_mmu_eg_db

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
    run_go_enrichment.R \\
        ${unfiltered_deg_flag} \\
        ${pvalue_threshold} \\
        ${logfc_threshold_flag} \\
        ${species_name} \\
        ${wgcna_list_flag} \\
        ${masigpro_list_flag} \\
        ${genes_diu_flag} \\
        ${genes_deu_flag} \\
        ${ntop_processes}
    """

    stub:
    """
    mkdir GO_results
    """
}
