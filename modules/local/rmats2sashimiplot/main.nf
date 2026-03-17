process RMATS2SASHIMIPLOT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e1/e11f3cff1d2216098d460c978eaeff3cfbc1847fc0c39e57eaf27966442c97b6/data'
        : 'community.wave.seqera.io/library/rmats2sashimiplot:3.0.0--efea77000abc79df'}"

    input:
    tuple val(meta), path(rmats_sigevent), path(bam_list1), path(bam_list2)
    tuple val(meta2), path(bam_files)

    output:
    tuple val(meta), path("${prefix}"), emit: sashimi_output
    tuple val("${task.process}"), val('rmats2sashimiplot'), eval("pip show rmats2sashimiplot | grep 'Version:' | awk '{print \$2}' | sed 's/Version: //'"), topic: versions, emit: versions_rmats2sashimiplot

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "sashimi_${meta.id}_out"
    if (!(meta.id in ["SE", "MXE", "A3SS", "A5SS", "RI"])) {
        error("Unknown event type: ${prefix}")
    }
    def condition_flag = bam_list1.baseName.toLowerCase() =~ /control/ ? "--l1 ${bam_list2.baseName} --l2 ${bam_list1.baseName}" : "--l1 ${bam_list1.baseName} --l2 ${bam_list2.baseName}"
    def bam_list_flag = bam_list1.baseName.toLowerCase() =~ /control/ ? "--b1 ${bam_list2} --b2 ${bam_list1}" : "--b1 ${bam_list1} --b2 ${bam_list2}"
    """
    rmats2sashimiplot \\
        -o ${prefix} \\
        ${condition_flag} \\
        --event-type ${meta.id} \\
        -e ${rmats_sigevent} \\
        ${bam_list_flag}
    """

    stub:
    prefix = task.ext.prefix ?: "sashimi_${meta.id}_out"
    """
    mkdir ${prefix}
    """
}
