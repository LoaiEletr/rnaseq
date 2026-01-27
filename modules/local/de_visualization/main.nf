#!/usr/bin/env nextflow

process DE_VISUALIZATION {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/72/72e333e1f91d258ebeb0c4c2f62d42037672bc6404138cb4ae9039dbc31c1d55/data'
        : 'community.wave.seqera.io/library/bioconductor-enhancedvolcano_r-base_r-pheatmap_r-rcolorbrewer:75c288f5aae6a1aa'}"

    input:
    path deg_file
    path samplesheet_csv
    path cluster_file
    val pval_threshold
    val lfc_threshold
    path vst_file

    output:
    path "de_visualization/volcano_sig_genes.pdf", optional: true, emit: volcano_plot
    path "de_visualization/heatmap/**", optional: true, emit: heatmaps
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def cluster_file_flag = cluster_file ?: null
    def deg_file_flag = deg_file ?: null
    """
    run_de_visualization.R \\
        ${deg_file_flag} \\
        ${samplesheet_csv} \\
        ${cluster_file_flag} \\
        ${pval_threshold} \\
        ${lfc_threshold} \\
        ${vst_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-enhancedvolcano: \$(Rscript -e "library(EnhancedVolcano); cat(as.character(packageVersion('EnhancedVolcano')))")
        r-pheatmap: \$(Rscript -e "library(pheatmap); cat(as.character(packageVersion('pheatmap')))")
        r-rcolorbrewer: \$(Rscript -e "library(RColorBrewer); cat(as.character(packageVersion('RColorBrewer')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p de_visualization/heatmap
    touch de_visualization/volcano_sig_genes.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-enhancedvolcano: \$(Rscript -e "library(EnhancedVolcano); cat(as.character(packageVersion('EnhancedVolcano')))")
        r-pheatmap: \$(Rscript -e "library(pheatmap); cat(as.character(packageVersion('pheatmap')))")
        r-rcolorbrewer: \$(Rscript -e "library(RColorBrewer); cat(as.character(packageVersion('RColorBrewer')))")
    END_VERSIONS
    """
}
