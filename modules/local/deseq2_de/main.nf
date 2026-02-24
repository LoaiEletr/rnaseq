process DESEQ2_DE {
    tag "${counts_rds}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7d/7d949caad176515ec71546b357488b0dce8cb84bea1c55ed71ea28c0d58c9383/data'
        : 'community.wave.seqera.io/library/bioconductor-deseq2_r-base:648c382534ca80b5'}"

    input:
    path counts_rds
    path samplesheet_csv
    val pvalue_threshold
    val logfc_threshold
    val top_genes

    output:
    path "unfiltered_counts.rds", emit: unfiltered_log2counts
    path "filtered_counts.rds", emit: filtered_log2counts
    path "normalized_counts.rds", emit: normalized_log2counts
    path "deseq2_de_genes.rds", emit: diffgenes_matrix
    path "vst_transformed.rds", emit: vst_object
    path "de_genes_full_table.rds", emit: diffgenes_table
    path "unfiltered_deg_results.rds", emit: unfiltered_deg
    path "top_${top_genes}_genes.rds", emit: topgenes_rds
    path "top_${top_genes}_genes.csv", emit: topgenes_csv
    tuple val("${task.process}"), val('r-base'), eval("Rscript -e 'cat(as.character(getRversion()))'"), topic: versions, emit: versions_rbase
    tuple val("${task.process}"), val('bioconductor-deseq2'), eval("Rscript -e \"library(DESeq2); cat(as.character(packageVersion('DESeq2')))\""), topic: versions, emit: versions_deseq2

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_deseq2_de.R \\
        ${counts_rds} \\
        ${samplesheet_csv} \\
        ${pvalue_threshold} \\
        ${logfc_threshold} \\
        ${top_genes}
    """

    stub:
    """
    touch unfiltered_counts.rds
    touch filtered_counts.rds
    touch normalized_counts.rds
    touch deseq2_de_genes.rds
    touch vst_transformed.rds
    touch de_genes_full_table.rds
    touch unfiltered_deg_results.rds
    touch top_${top_genes}_genes.rds
    touch top_${top_genes}_genes.csv
    """
}
