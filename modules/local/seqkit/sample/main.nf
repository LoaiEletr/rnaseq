process SEQKIT_SAMPLE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/85/85b40b925e4d4a62f9b833bbb0646d7ea6cf53d8a875e3055f90da757d7ccd27/data'
        : 'community.wave.seqera.io/library/seqkit:2.12.0--5e60eb93d3a59212'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_1M.fastq.gz"), emit: reads
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/seqkit v//'"), topic: versions, emit: versions_seqkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        seqkit \\
            sample \\
            -s100 \\
            ${reads} \\
            -n 1000000 \\
            ${args} \\
            -o ${prefix}_1M.fastq.gz
        """
    }
    else {
        """
        seqkit \\
            sample \\
            -s100 \\
            ${reads[0]} \\
            -n 1000000 \\
            ${args} \\
            -o ${prefix}_1_1M.fastq.gz

        seqkit \\
            sample \\
            -s100 \\
            ${reads[1]} \\
            -n 1000000 \\
            -j ${task.cpus} \\
            ${args} \\
            -o ${prefix}_2_1M.fastq.gz
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        echo "" | gzip > ${prefix}_1M.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$( seqkit version | sed 's/seqkit v//' )
        END_VERSIONS
        """
    }
    else {
        """
        echo "" | gzip > ${prefix}_1_1M.fastq.gz
        echo "" | gzip > ${prefix}_2_1M.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$( seqkit version | sed 's/seqkit v//' )
        END_VERSIONS
        """
    }
}
