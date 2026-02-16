process BCFTOOLS_INDEX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5b/5bfed643f44786895efe9c8894620b3a41c597f07b7aa7c51d22cd3f66411411/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23--c4af0e3be58276b3'}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf.gz.csi"), emit: csi
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bcftools \\
        index \\
        ${args} \\
        --threads ${task.cpus} \\
        ${vcf}
    """

    stub:
    """
    touch ${vcf}.csi
    """
}
