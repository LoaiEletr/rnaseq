process TX2GENE {
    tag "${species_name}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f5/f58661c83be1d4605d0902036012cf776c0744f8fdc1ec3c5d8dcaed038b0553/data'
        : 'community.wave.seqera.io/library/bioconductor-biomart_r-base:f9e5215e8454b7b6'}"

    input:
    val species_name

    output:
    path "*.tsv", emit: tsv
    tuple val("${task.process}"), val('r-base'), eval("Rscript -e 'cat(as.character(getRversion()))'"), topic: versions, emit: versions_rbase
    tuple val("${task.process}"), val('bioconductor-biomart'), eval("Rscript -e \"cat(as.character(packageVersion('biomaRt')))\""), topic: versions, emit: versions_biomart

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    generate_tx2gene.R ${species_name}
    """

    stub:
    """
    touch tx2gene.tsv
    """
}
