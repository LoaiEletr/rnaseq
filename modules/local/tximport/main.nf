process TXIMPORT {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b8/b8617b62baac3163f0458b3f72ede655a6376c5b713f7d7b7507c6dbe7209f9e/data'
        : 'community.wave.seqera.io/library/bioconductor-rhdf5_bioconductor-tximport_r-base_r-jsonlite:0454e92595745a41'}"

    input:
    path quant_dirs
    val quant_type
    path tx2gene_table

    output:
    path "*.rds", emit: rds
    tuple val("${task.process}"), val('r-base'), eval("Rscript -e 'cat(as.character(getRversion()))'"), topic: versions, emit: versions_rbase
    tuple val("${task.process}"), val('bioconductor-tximport'), eval("Rscript -e \"cat(as.character(packageVersion('tximport')))\""), topic: versions, emit: versions_tximport
    tuple val("${task.process}"), val('bioconductor-rhdf5'), eval("Rscript -e \"cat(as.character(packageVersion('rhdf5')))\""), topic: versions, emit: versions_rhdf5

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_tximport.R ${quant_dirs} ${quant_type} ${tx2gene_table}
    """

    stub:
    """
    touch tximport_results.rds
    """
}
