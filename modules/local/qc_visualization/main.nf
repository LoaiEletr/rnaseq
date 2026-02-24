process QC_VISUALIZATION {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c3/c3595ff0bc3d232084c47e0ab7b432ddccac7ba9213bce838fb816425ba646a7/data'
        : 'community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-edger_r-base_r-cowplot_pruned:ccc132c22a3d4a51'}"

    input:
    path voom_transformed
    path samplesheet_csv
    path unfiltered_log2counts
    path filtered_logcounts
    path normalized_logcounts
    path vst_object

    output:
    path "qc_plots", emit: qc_plots
    tuple val("${task.process}"), val('r-base'), eval("Rscript -e 'cat(as.character(getRversion()))'"), topic: versions, emit: versions_rbase
    tuple val("${task.process}"), val('bioconductor-edger'), eval("Rscript -e \"cat(as.character(packageVersion('edgeR')))\""), topic: versions, emit: versions_edger
    tuple val("${task.process}"), val('r-ggplot2'), eval("Rscript -e \"cat(as.character(packageVersion('ggplot2')))\""), topic: versions, emit: versions_ggplot2
    tuple val("${task.process}"), val('bioconductor-deseq2'), eval("Rscript -e \"cat(as.character(packageVersion('DESeq2')))\""), topic: versions, emit: versions_deseq2
    tuple val("${task.process}"), val('r-tibble'), eval("Rscript -e \"cat(as.character(packageVersion('tibble')))\""), topic: versions, emit: versions_tibble
    tuple val("${task.process}"), val('r-tidyr'), eval("Rscript -e \"cat(as.character(packageVersion('tidyr')))\""), topic: versions, emit: versions_tidyr
    tuple val("${task.process}"), val('r-dplyr'), eval("Rscript -e \"cat(as.character(packageVersion('dplyr')))\""), topic: versions, emit: versions_dplyr
    tuple val("${task.process}"), val('r-cowplot'), eval("Rscript -e \"cat(as.character(packageVersion('cowplot')))\""), topic: versions, emit: versions_cowplot
    tuple val("${task.process}"), val('r-mass'), eval("Rscript -e \"cat(as.character(packageVersion('MASS')))\""), topic: versions, emit: versions_mass

    when:
    task.ext.when == null || task.ext.when

    script:
    def voom_transformed_flag = voom_transformed ?: null
    def vst_object_flag = vst_object ?: normalized_logcounts
    """
    run_qc_visualization.R \\
        ${voom_transformed_flag} \\
        ${samplesheet_csv} \\
        ${unfiltered_log2counts} \\
        ${filtered_logcounts} \\
        ${normalized_logcounts} \\
        ${vst_object_flag}
    """

    stub:
    """
    mkdir -p qc_plots
    touch qc_plots/mds_plot.pdf
    touch qc_plots/pca_plot.pdf
    touch qc_plots/violin_plots.pdf
    """
}
