process RSEQC_GENEBODYCOVERAGE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1d/1db6950626c14db8a5b5a80089c259774a693fc3a5946d1bf169d19b11f7bccb/data'
        : 'community.wave.seqera.io/library/rseqc_r-base:a2f5852c8ab06f43'}"

    input:
    tuple val(meta), path(bam_files)
    tuple val(meta2), path(bai_files)
    tuple val(meta3), path(bed)

    output:
    tuple val(meta), path("*.geneBodyCoverage.r"), emit: rscript
    tuple val(meta), path("*.geneBodyCoverage.curves.pdf"), emit: curves
    tuple val(meta), path("*.geneBodyCoverage.heatMap.pdf"), optional: true, emit: heatmap
    tuple val(meta), path("*.geneBodyCoverage.txt"), emit: txt
    tuple val(meta), path("log.txt"), emit: log
    tuple val("${task.process}"), val('rseqc'), eval("geneBody_coverage.py --version | sed 's/.*.py //'"), topic: versions, emit: versions_rseqc

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "genebodycoverage"
    """
    geneBody_coverage.py \\
        -r ${bed} \\
        -i ${bam_files.join(',')} \\
        -o ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "genebodycoverage"
    """
    touch ${prefix}.geneBodyCoverage.r
    touch ${prefix}.geneBodyCoverage.curves.pdf
    touch ${prefix}.geneBodyCoverage.heatMap.pdf
    touch ${prefix}.geneBodyCoverage.txt
    touch log.txt
    """
}
