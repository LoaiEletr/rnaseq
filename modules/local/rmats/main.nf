#!/usr/bin/env nextflow

process RMATS {
    tag "${mode}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/rmats:4.3.0--py39hbadf43b_5'
        : 'biocontainers/rmats:4.3.0--py39hbadf43b_5'}"

    input:
    tuple val(meta), path(bam_files)
    path bam_list
    path gtf
    val mode
    val avg_length
    val statoff
    val novelss
    val allowclipping
    val individualcounts
    path rmats_prep_tmp_dir

    output:
    path ("rmats_post_output"), emit: post_output, optional: true
    path ("rmats_prep_tmp"), emit: prep_tmp, optional: true
    path ("rmats_prep_output"), emit: prep_output, optional: true
    path ("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (!(mode in ["prep", "post"])) {
        error("Unknown mode: ${mode}")
    }
    def prefix = mode == 'post' ? 'rmats_post_output' : 'rmats_prep_output'
    def bam_list_flag = bam_list[0].baseName.toLowerCase() =~ /control/ ? "--b1 ${bam_list[1]} --b2 ${bam_list[0]}" : "--b1 ${bam_list[0]} --b2 ${bam_list[1]}"
    def rmats_output = "--od ${prefix}"
    def rmats_tmp = mode == 'post' && "${rmats_prep_tmp_dir}" ? "--tmp ${rmats_prep_tmp_dir}" : '--tmp rmats_prep_tmp'
    def statoff_flag = statoff ? '--statoff' : ''
    def novelSS_flag = novelss ? '--novelSS' : ''
    def allowClipping_flag = allowclipping ? '--allow-clipping' : ''
    def individualCounts_flag = individualcounts ? '--individual-counts' : ''
    read_type = meta.single_end ? "-t single" : "-t paired"
    strand_flag = ""
    if (meta.lib_type == "IU" || meta.lib_type == "U") {
        strand_flag = "--libType fr-unstranded"
    }
    else if (meta.lib_type in ["SR", "ISR", "MSR", "OSR", "reverse"]) {
        strand_flag = "--libType fr-firststrand"
    }
    else if (meta.lib_type in ["SF", "ISF", "MSF", "OSF", "forward"]) {
        strand_flag = "--libType fr-secondstrand"
    }
    else {
        error("Unknown strandness type: ${meta.lib_type}")
    }
    """
    rmats.py \\
        --gtf ${gtf} \\
        ${bam_list_flag} \\
        ${rmats_output} \\
        ${rmats_tmp} \\
        ${read_type} \\
        ${strand_flag} \\
        --readLength ${avg_length} \\
        --variable-read-length \\
        --nthread ${task.cpus} \\
        --task ${mode} \\
        ${statoff_flag} \\
        ${novelSS_flag} \\
        ${allowClipping_flag} \\
        ${individualCounts_flag} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rmats: \$(rmats.py --version | sed 's/v//')
    END_VERSIONS
    """

    stub:
    def prefix = mode == 'post' ? 'rmats_post_output' : 'rmats_prep_output'
    """
    # Create output directory
    mkdir -p ${prefix}

    # Create all rMATS output files
    touch ${prefix}/A3SS.MATS.JCEC.txt
    touch ${prefix}/A3SS.MATS.JC.txt
    touch ${prefix}/A5SS.MATS.JCEC.txt
    touch ${prefix}/A5SS.MATS.JC.txt
    touch ${prefix}/MXE.MATS.JCEC.txt
    touch ${prefix}/MXE.MATS.JC.txt
    touch ${prefix}/RI.MATS.JCEC.txt
    touch ${prefix}/RI.MATS.JC.txt
    touch ${prefix}/SE.MATS.JCEC.txt
    touch ${prefix}/SE.MATS.JC.txt

    touch ${prefix}/fromGTF.A3SS.txt
    touch ${prefix}/fromGTF.A5SS.txt
    touch ${prefix}/fromGTF.MXE.txt
    touch ${prefix}/fromGTF.RI.txt
    touch ${prefix}/fromGTF.SE.txt

    touch ${prefix}/fromGTF.novelJunction.A3SS.txt
    touch ${prefix}/fromGTF.novelJunction.A5SS.txt
    touch ${prefix}/fromGTF.novelJunction.MXE.txt
    touch ${prefix}/fromGTF.novelJunction.RI.txt
    touch ${prefix}/fromGTF.novelJunction.SE.txt

    touch ${prefix}/fromGTF.novelSpliceSite.A3SS.txt
    touch ${prefix}/fromGTF.novelSpliceSite.A5SS.txt
    touch ${prefix}/fromGTF.novelSpliceSite.MXE.txt
    touch ${prefix}/fromGTF.novelSpliceSite.RI.txt
    touch ${prefix}/fromGTF.novelSpliceSite.SE.txt

    touch ${prefix}/JCEC.raw.input.A3SS.txt
    touch ${prefix}/JCEC.raw.input.A5SS.txt
    touch ${prefix}/JCEC.raw.input.MXE.txt
    touch ${prefix}/JCEC.raw.input.RI.txt
    touch ${prefix}/JCEC.raw.input.SE.txt

    touch ${prefix}/JC.raw.input.A3SS.txt
    touch ${prefix}/JC.raw.input.A5SS.txt
    touch ${prefix}/JC.raw.input.MXE.txt
    touch ${prefix}/JC.raw.input.RI.txt
    touch ${prefix}/JC.raw.input.SE.txt

    touch ${prefix}/summary.txt

    # Create tmp directory with its files
    mkdir -p ${prefix}/tmp/JC_A3SS
    touch ${prefix}/tmp/JC_A3SS/rMATS_result_FDR.txt
    touch ${prefix}/tmp/JC_A3SS/rMATS_result_ID.txt
    touch ${prefix}/tmp/JC_A3SS/rMATS_result_I-L.txt
    touch ${prefix}/tmp/JC_A3SS/rMATS_result_INP.txt
    touch ${prefix}/tmp/JC_A3SS/rMATS_result_P-V.txt
    touch ${prefix}/tmp/JC_A3SS/rMATS_result_.txt

    mkdir -p ${prefix}/tmp/JC_A5SS
    touch ${prefix}/tmp/JC_A5SS/rMATS_result_FDR.txt
    touch ${prefix}/tmp/JC_A5SS/rMATS_result_ID.txt
    touch ${prefix}/tmp/JC_A5SS/rMATS_result_I-L.txt
    touch ${prefix}/tmp/JC_A5SS/rMATS_result_INP.txt
    touch ${prefix}/tmp/JC_A5SS/rMATS_result_P-V.txt
    touch ${prefix}/tmp/JC_A5SS/rMATS_result_.txt

    mkdir -p ${prefix}/tmp/JC_MXE
    touch ${prefix}/tmp/JC_MXE/rMATS_result_FDR.txt
    touch ${prefix}/tmp/JC_MXE/rMATS_result_ID.txt
    touch ${prefix}/tmp/JC_MXE/rMATS_result_I-L.txt
    touch ${prefix}/tmp/JC_MXE/rMATS_result_INP.txt
    touch ${prefix}/tmp/JC_MXE/rMATS_result_P-V.txt
    touch ${prefix}/tmp/JC_MXE/rMATS_result_.txt

    mkdir -p ${prefix}/tmp/JC_RI
    touch ${prefix}/tmp/JC_RI/rMATS_result_FDR.txt
    touch ${prefix}/tmp/JC_RI/rMATS_result_ID.txt
    touch ${prefix}/tmp/JC_RI/rMATS_result_I-L.txt
    touch ${prefix}/tmp/JC_RI/rMATS_result_INP.txt
    touch ${prefix}/tmp/JC_RI/rMATS_result_P-V.txt
    touch ${prefix}/tmp/JC_RI/rMATS_result_.txt

    mkdir -p ${prefix}/tmp/JC_SE
    touch ${prefix}/tmp/JC_SE/rMATS_result_FDR.txt
    touch ${prefix}/tmp/JC_SE/rMATS_result_ID.txt
    touch ${prefix}/tmp/JC_SE/rMATS_result_I-L.txt
    touch ${prefix}/tmp/JC_SE/rMATS_result_INP.txt
    touch ${prefix}/tmp/JC_SE/rMATS_result_P-V.txt
    touch ${prefix}/tmp/JC_SE/rMATS_result_.txt

    mkdir -p ${prefix}/tmp/JCEC_A3SS
    touch ${prefix}/tmp/JCEC_A3SS/rMATS_result_FDR.txt
    touch ${prefix}/tmp/JCEC_A3SS/rMATS_result_ID.txt
    touch ${prefix}/tmp/JCEC_A3SS/rMATS_result_I-L.txt
    touch ${prefix}/tmp/JCEC_A3SS/rMATS_result_INP.txt
    touch ${prefix}/tmp/JCEC_A3SS/rMATS_result_P-V.txt
    touch ${prefix}/tmp/JCEC_A3SS/rMATS_result_.txt

    mkdir -p ${prefix}/tmp/JCEC_A5SS
    touch ${prefix}/tmp/JCEC_A5SS/rMATS_result_FDR.txt
    touch ${prefix}/tmp/JCEC_A5SS/rMATS_result_ID.txt
    touch ${prefix}/tmp/JCEC_A5SS/rMATS_result_I-L.txt
    touch ${prefix}/tmp/JCEC_A5SS/rMATS_result_INP.txt
    touch ${prefix}/tmp/JCEC_A5SS/rMATS_result_P-V.txt
    touch ${prefix}/tmp/JCEC_A5SS/rMATS_result_.txt

    mkdir -p ${prefix}/tmp/JCEC_MXE
    touch ${prefix}/tmp/JCEC_MXE/rMATS_result_FDR.txt
    touch ${prefix}/tmp/JCEC_MXE/rMATS_result_ID.txt
    touch ${prefix}/tmp/JCEC_MXE/rMATS_result_I-L.txt
    touch ${prefix}/tmp/JCEC_MXE/rMATS_result_INP.txt
    touch ${prefix}/tmp/JCEC_MXE/rMATS_result_P-V.txt
    touch ${prefix}/tmp/JCEC_MXE/rMATS_result_.txt

    mkdir -p ${prefix}/tmp/JCEC_RI
    touch ${prefix}/tmp/JCEC_RI/rMATS_result_FDR.txt
    touch ${prefix}/tmp/JCEC_RI/rMATS_result_ID.txt
    touch ${prefix}/tmp/JCEC_RI/rMATS_result_I-L.txt
    touch ${prefix}/tmp/JCEC_RI/rMATS_result_INP.txt
    touch ${prefix}/tmp/JCEC_RI/rMATS_result_P-V.txt
    touch ${prefix}/tmp/JCEC_RI/rMATS_result_.txt

    mkdir -p ${prefix}/tmp/JCEC_SE
    touch ${prefix}/tmp/JCEC_SE/rMATS_result_FDR.txt
    touch ${prefix}/tmp/JCEC_SE/rMATS_result_ID.txt
    touch ${prefix}/tmp/JCEC_SE/rMATS_result_I-L.txt
    touch ${prefix}/tmp/JCEC_SE/rMATS_result_INP.txt
    touch ${prefix}/tmp/JCEC_SE/rMATS_result_P-V.txt
    touch ${prefix}/tmp/JCEC_SE/rMATS_result_.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rmats: \$(rmats.py --version | sed 's/v//')
    END_VERSIONS
    """
}
