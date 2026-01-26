#!/usr/bin/env nextflow

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
    path "versions.yml", emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-edger: \$(Rscript -e "library(edgeR); cat(as.character(packageVersion('edgeR')))")
        r-ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
        r-tibble: \$(Rscript -e "library(tibble); cat(as.character(packageVersion('tibble')))")
        r-tidyr: \$(Rscript -e "library(tidyr); cat(as.character(packageVersion('tidyr')))")
        r-dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        r-cowplot: \$(Rscript -e "library(cowplot); cat(as.character(packageVersion('cowplot')))")
        r-mass: \$(Rscript -e "library(MASS); cat(as.character(packageVersion('MASS')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p qc_plots
    touch qc_plots/mds_plot.pdf
    touch qc_plots/pca_plot.pdf
    touch qc_plots/violin_plots.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-edger: \$(Rscript -e "library(edgeR); cat(as.character(packageVersion('edgeR')))")
        r-ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')))")
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
        r-tibble: \$(Rscript -e "library(tibble); cat(as.character(packageVersion('tibble')))")
        r-tidyr: \$(Rscript -e "library(tidyr); cat(as.character(packageVersion('tidyr')))")
        r-dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
        r-cowplot: \$(Rscript -e "library(cowplot); cat(as.character(packageVersion('cowplot')))")
        r-mass: \$(Rscript -e "library(MASS); cat(as.character(packageVersion('MASS')))")
    END_VERSIONS
    """
}
