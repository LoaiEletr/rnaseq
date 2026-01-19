#!/usr/bi/env nextflow

process RSEQC_GENEBODYCOVERAGE {
    tag "${bam_files.size()} samples"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-c80ae8d0fe5685926c9bc673e400ff09a71844fd:29c8e89bc12d33b39e760c5ca3b1cfa087927580-0'
        : 'biocontainers/mulled-v2-c80ae8d0fe5685926c9bc673e400ff09a71844fd:e01414b01cd5729e641b84adaf1ce4fd6181bcb8-2'}"

    input:
    path bam_files
    path bai_files
    path bed

    output:
    path "*.geneBodyCoverage.r", emit: rscript
    path "*.geneBodyCoverage.curves.pdf", emit: pdf
    path "*.geneBodyCoverage.txt", emit: txt
    path "log.txt", emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "genebody_coverage"
    """
    geneBody_coverage.py \\
        -r ${bed} \\
        -i ${bam_files.join(',')} \\
        -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$( geneBody_coverage.py --version | sed 's/.*.py //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "genebody_coverage"
    """
    touch ${prefix}.geneBodyCoverage.r
    touch ${prefix}.geneBodyCoverage.curves.pdf
    touch ${prefix}.geneBodyCoverage.txt
    touch log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$( geneBody_coverage.py --version | sed 's/.*.py //' )
    END_VERSIONS
    """
}
