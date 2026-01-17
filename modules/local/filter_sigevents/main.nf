#!/usr/bin/env nextflow

process FILTER_SIGEVENTS {
    tag "${event_type}"
    label 'process_low'

    input:
    path rmats_post_output
    each event_type
    val fdr_cutoff
    val delta_psi

    output:
    tuple val(event_type), path("${event_type}_filtered.txt"), emit: sig_rmats
    path ("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${event_type}_filtered.txt"
    """
    cat ${rmats_post_output}/${event_type}.MATS.JC.txt | awk -v fdr=${fdr_cutoff} -v dpsi=${delta_psi} '\\
    NR==1 { print; next } \\
    { \\
        split(\$13,a,","); split(\$14,b,","); \\
        sum1=0; for(i in a) sum1+=a[i]; for(i in b) sum1+=b[i]; avg1=sum1/(length(a)+length(b)); \\
        split(\$15,c,","); split(\$16,d,","); \\
        sum2=0; for(i in c) sum2+=c[i]; for(i in d) sum2+=d[i]; avg2=sum2/(length(c)+length(d)); \\
        split(\$21,p,","); sump=0; countp=0; for(i in p){ if(p[i]!="NA"){ sump+=p[i]; countp++ } }; avgp=sump/countp; \\
        split(\$22,q,","); sumq=0; countq=0; for(i in q){ if(q[i]!="NA"){ sumq+=q[i]; countq++ } }; avgq=sumq/countq; \\
        if (\$20 <= fdr && (\$23 >= dpsi || \$23 <= -dpsi) && avg1 >= 10 && avg2 >= 10 && avgp >= 0.05 && avgp <= 0.95 && avgq >= 0.05 && avgq <= 0.95) \\
            print; \\
    }' > ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | awk 'NR==1 {print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${event_type}_filtered.txt"
    """
    touch ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | awk 'NR==1 {print \$2}')
    END_VERSIONS
    """
}
