process KALLISTO_QUANT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/70/70a7a1b2f5dc1ca0e801e958ab14ea029f662f1ead026ba9cdff59f99995f19c/data'
        : 'community.wave.seqera.io/library/kallisto:0.51.1--b63691b6841c7a52'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(gtf)
    val fragment_length
    val fragment_length_sd

    output:
    tuple val(meta), path("${prefix}"), emit: quant_dir
    tuple val(meta), path("${prefix}/*.log"), emit: log
    tuple val("${task.process}"), val('kallisto'), eval("kallisto version | sed 's/.*ion //'"), topic: versions, emit: versions_kallisto

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def fastq_in = meta.single_end ? "--single -l ${fragment_length} -s ${fragment_length_sd} ${reads}" : "${reads[0]} ${reads[1]}"
    def strand_flag = ""

    if (meta.lib_type in ["IU", "U", "MU", "OU"]) {
        strand_flag = ""
    }
    else if (meta.lib_type in ["ISR", "SR", "MSR", "OSR", "reverse"]) {
        strand_flag = "--rf-stranded"
    }
    else if (meta.lib_type in ["ISF", "SF", "MSF", "OSF", "forward"]) {
        strand_flag = "--fr-stranded"
    }
    else {
        error("Unknown meta.lib_type type: ${meta.lib_type}")
    }
    """
    kallisto \\
        quant \\
        ${args} \\
        --gtf ${gtf} \\
        ${strand_flag} \\
        -t ${task.cpus} \\
        -i ${index} \\
        -o ${prefix} \\
        ${fastq_in} \\
        2>| >( tee kallisto_${prefix}_quant.log >&2 )

    mv kallisto_${prefix}_quant.log ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}

    touch kallisto_${prefix}_quant.log
    touch ${prefix}/abundance.h5
    touch ${prefix}/abundance.tsv
    touch ${prefix}/run_info.json

    mv kallisto_${prefix}_quant.log ${prefix}
    """
}
