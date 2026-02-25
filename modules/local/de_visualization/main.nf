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
    tuple val("${task.process}"), val('r-base'), eval("Rscript -e 'cat(as.character(getRversion()))'"), topic: versions, emit: versions_rbase
    tuple val("${task.process}"), val('bioconductor-enhancedvolcano'), eval("Rscript -e \"library(EnhancedVolcano); cat(as.character(packageVersion('EnhancedVolcano')))\""), topic: versions, emit: versions_enhancedvolcano
    tuple val("${task.process}"), val('r-pheatmap'), eval("Rscript -e \"library(pheatmap); cat(as.character(packageVersion('pheatmap')))\""), topic: versions, emit: versions_pheatmap
    tuple val("${task.process}"), val('r-rcolorbrewer'), eval("Rscript -e \"library(RColorBrewer); cat(as.character(packageVersion('RColorBrewer')))\""), topic: versions, emit: versions_rcolorbrewer

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
    """

    stub:
    """
    mkdir -p de_visualization/heatmap
    touch de_visualization/volcano_sig_genes.pdf
    """
}
