process MSIGDB_GSEA {
    tag "${unfiltered_deg}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9e/9e4333b647c7183896fe4cd6054c7733a9a5b93771d50acda8a5169436591e02/data'
        : 'community.wave.seqera.io/library/bioconductor-clusterprofiler_bioconductor-enrichplot_r-base_r-dplyr_pruned:26af7e9ffc04c073'}"

    input:
    path unfiltered_deg
    val pvalue_threshold
    val species_name
    val msigdb_categories
    val nes_threshold
    val padj_gsea
    val ntop_processes
    val rank_method

    output:
    path "GSEA_results", emit: gsea_results, optional: true
    tuple val("${task.process}"), val('r-base'), eval("Rscript -e 'cat(as.character(getRversion()))'"), topic: versions, emit: versions_rbase
    tuple val("${task.process}"), val('bioconductor-clusterprofiler'), eval("Rscript -e \"cat(as.character(packageVersion('clusterProfiler')))\""), topic: versions, emit: versions_clusterprofiler
    tuple val("${task.process}"), val('bioconductor-enrichplot'), eval("Rscript -e \"cat(as.character(packageVersion('enrichplot')))\""), topic: versions, emit: versions_enrichplot
    tuple val("${task.process}"), val('r-ggplot2'), eval("Rscript -e \"cat(as.character(packageVersion('ggplot2')))\""), topic: versions, emit: versions_ggplot2
    tuple val("${task.process}"), val('r-dplyr'), eval("Rscript -e \"cat(as.character(packageVersion('dplyr')))\""), topic: versions, emit: versions_dplyr
    tuple val("${task.process}"), val('r-tibble'), eval("Rscript -e \"cat(as.character(packageVersion('tibble')))\""), topic: versions, emit: versions_tibble
    tuple val("${task.process}"), val('r-msigdbr'), eval("Rscript -e \"cat(as.character(packageVersion('msigdbr')))\""), topic: versions, emit: versions_msigdbr

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_msigdb_gsea.R \\
        ${unfiltered_deg} \\
        ${pvalue_threshold} \\
        ${species_name} \\
        ${msigdb_categories} \\
        ${nes_threshold} \\
        ${padj_gsea} \\
        ${ntop_processes} \\
        ${rank_method}
    """

    stub:
    """
    mkdir GSEA_results
    """
}
