process BCFTOOLS_ISEC {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5b/5bfed643f44786895efe9c8894620b3a41c597f07b7aa7c51d22cd3f66411411/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23--c4af0e3be58276b3'}"

    input:
    tuple val(meta), path(vcfs), path(tbis)
    tuple val(meta2), path(file_list)
    tuple val(met3), path(regions_file)

    output:
    tuple val(meta), path("${prefix}", type: "dir"), emit: results
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def regions_file_flag = regions_file ? "-R ${regions_file}" : ''
    def vcf_files_flag = file_list ? "-l ${file_list}" : "${vcfs}"
    def vcf0_base = vcfs[0].name.tokenize('.')[0]
    def vcf1_base = vcfs[1].name.tokenize('.')[0]
    def rename_command = file_list
        ? ''
        : """
    # 1. Update the internal text of the README so the names match your new files
    if [ -f ${prefix}/README.txt ]; then
        sed -i "s|0000.vcf.gz|${vcf0_base}_unique.vcf.gz|g" ${prefix}/README.txt
        sed -i "s|0001.vcf.gz|${vcf1_base}_unique.vcf.gz|g" ${prefix}/README.txt
        sed -i "s|0002.vcf.gz|${vcf0_base}_common.vcf.gz|g" ${prefix}/README.txt
        sed -i "s|0003.vcf.gz|${vcf1_base}_common.vcf.gz|g" ${prefix}/README.txt
    fi

    # 2. Rename the actual VCF files and their TBI indices
    mv ${prefix}/0000.vcf.gz ${prefix}/${vcf0_base}_unique.vcf.gz
    mv ${prefix}/0000.vcf.gz.tbi ${prefix}/${vcf0_base}_unique.vcf.gz.tbi

    mv ${prefix}/0001.vcf.gz ${prefix}/${vcf1_base}_unique.vcf.gz
    mv ${prefix}/0001.vcf.gz.tbi ${prefix}/${vcf1_base}_unique.vcf.gz.tbi

    mv ${prefix}/0002.vcf.gz ${prefix}/${vcf0_base}_common.vcf.gz
    mv ${prefix}/0002.vcf.gz.tbi ${prefix}/${vcf0_base}_common.vcf.gz.tbi

    mv ${prefix}/0003.vcf.gz ${prefix}/${vcf1_base}_common.vcf.gz
    mv ${prefix}/0003.vcf.gz.tbi ${prefix}/${vcf1_base}_common.vcf.gz.tbi
    """.stripIndent()

    """
    bcftools isec  \\
        ${args} \\
        ${regions_file_flag} \\
        -p ${prefix} \\
        ${vcf_files_flag}

    ${rename_command}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/README.txt
    touch ${prefix}/sites.txt
    echo "" | gzip > ${prefix}/0000.vcf.gz
    touch ${prefix}/0000.vcf.gz.tbi
    echo "" | gzip > ${prefix}/0001.vcf.gz
    touch ${prefix}/0001.vcf.gz.tbi
    """
}
