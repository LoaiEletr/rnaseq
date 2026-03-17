//
// Subsample FASTQ reads, quantify with Salmon, and determine strandness
// Performs: SeqKit subsampling, Salmon quantification and strandness detection
//

include { SEQKIT_SAMPLE } from '../../../modules/local/seqkit/sample/main.nf'
include { SALMON_QUANT } from '../../../modules/local/salmon/quant/main.nf'

workflow FASTQ_SAMPLE_SEQKIT_SALMON {
    take:
    ch_reads // channel: [ val(meta), reads ]
    ch_gtf // channel: [ val(meta), gtf ]
    ch_salmon_index // channel: [ val(meta), index ]

    main:

    // 1. Subsample reads with SeqKit
    SEQKIT_SAMPLE(ch_reads)

    // 2. Quantify with Salmon and determine strandness
    SALMON_QUANT(
        SEQKIT_SAMPLE.out.reads,
        ch_salmon_index,
        ch_gtf,
    )

    emit:
    reads = SEQKIT_SAMPLE.out.reads // channel: [ val(meta), [ reads ] ]
    quant_dir = SALMON_QUANT.out.quant_dir // channel: [ val(meta), [ quant_dir ] ]
    strandness = SALMON_QUANT.out.strandness // channel: [ val(meta), [ strandness ] ]
}
