process ISOFORMSWITCHANALYZER {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3e/3e73e76384587e47c557b88bfce8212a821888d52b4647967e69895fe64f2083/data'
        : 'community.wave.seqera.io/library/bioconductor-isoformswitchanalyzer_bioconductor-rhdf5_r-base_r-ggplot2_r-mass:e41cc44dc5deb629'}"

    input:
    path quant_dirs
    path samplesheet
    val quant_type
    tuple val(meta), path(gtf)
    tuple val(meta2), path(transcript)
    val method
    val dif_cutoff
    val pval_threshold
    val ntop_isoforms

    output:
    path "isoform_output", emit: isoform_results, optional: true
    path "isoform_output/significant_DIU_gene_ids.rds", emit: significant_geneids, optional: true
    tuple val("${task.process}"), val('r-base'), eval("Rscript -e 'cat(as.character(getRversion()))'"), topic: versions, emit: versions_rbase
    tuple val("${task.process}"), val('bioconductor-isoformswitchanalyzer'), eval("Rscript -e \"cat(as.character(packageVersion('IsoformSwitchAnalyzeR')))\""), topic: versions, emit: versions_isoformswitchanalyzer
    tuple val("${task.process}"), val('bioconductor-rhdf5'), eval("Rscript -e \"cat(as.character(packageVersion('rhdf5')))\""), topic: versions, emit: versions_rhdf5
    tuple val("${task.process}"), val('r-mass'), eval("Rscript -e \"cat(as.character(packageVersion('MASS')))\""), topic: versions, emit: versions_mass
    tuple val("${task.process}"), val('r-ggplot2'), eval("Rscript -e \"cat(as.character(packageVersion('ggplot2')))\""), topic: versions, emit: versions_ggplot2

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_isoform_switch_analysis.R \\
        ${quant_dirs} \\
        ${samplesheet} \\
        ${quant_type} \\
        ${gtf} \\
        ${transcript} \\
        ${method} \\
        ${dif_cutoff} \\
        ${pval_threshold} \\
        ${ntop_isoforms}
    """

    stub:
    """
    mkdir isoform_output
    touch isoform_output/significant_DIU_gene_ids.rds
    """
}
