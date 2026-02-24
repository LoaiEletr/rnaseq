process LIMMA_DE {
    tag "${voom_object}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd1284a73c2a23f3cef38ec39f1695005012adab9c3ecc4ba56ed33a0b4daa14/data'
        : 'community.wave.seqera.io/library/bioconductor-limma_r-base:2579913b34b8b4da'}"

    input:
    path voom_object
    path samplesheet
    val pvalue_threshold
    val logfc_threshold
    val ntop_genes

    output:
    path "ebfit.rds", emit: ebfit_rds
    path "de_expression_values.rds", emit: de_expression_rds
    path "de_genes_full_table.rds", emit: de_genes_rds
    path "unfiltered_deg_results.rds", emit: unfiltered_deg_rds
    path "top_${ntop_genes}_genes.rds", emit: topgenes_rds
    path "top_${ntop_genes}_genes.csv", emit: topgenes_csv
    tuple val("${task.process}"), val('r-base'), eval("Rscript -e 'cat(as.character(getRversion()))'"), topic: versions, emit: versions_rbase
    tuple val("${task.process}"), val('bioconductor-limma'), eval("Rscript -e \"cat(as.character(packageVersion('limma')))\""), topic: versions, emit: versions_limma

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_limma_de.R ${voom_object} ${samplesheet} ${pvalue_threshold} ${logfc_threshold} ${ntop_genes}
    """

    stub:
    """
    touch ebfit.rds
    touch de_expression_values.rds
    touch de_genes_full_table.rds
    touch unfiltered_deg_results.rds
    touch top_${ntop_genes}_genes.rds
    touch top_${ntop_genes}_genes.csv
    """
}
