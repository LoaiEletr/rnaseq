//
// Subsample FASTQ reads, quantify with Salmon, and determine strandness
// Performs: SeqKit subsampling, Salmon indexing, Salmon quantification and strandness detection
//

include { SEQKIT_SAMPLE } from '../../../modules/local/seqkit/sample/main.nf'
include { SALMON_INDEX } from '../../../modules/local/salmon/index/main.nf'
include { SALMON_QUANT } from '../../../modules/local/salmon/quant/main.nf'

workflow FASTQ_SAMPLE_SEQKIT_SALMON {
    take:
    ch_reads // channel: [ val(meta), reads ]
    ch_transcriptome // channel: [ transcriptome ]
    ch_gtf // channel: [ gtf ]

    main:

    ch_versions = channel.empty()

    // 1. Subsample reads with SeqKit
    SEQKIT_SAMPLE(ch_reads)
    ch_versions = ch_versions.mix(SEQKIT_SAMPLE.out.versions.first())

    // 2. Build Salmon index
    SALMON_INDEX(ch_transcriptome)
    ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)

    // 3. Quantify with Salmon and determine strandness
    def libtype = 'A'
    // Automatic library type detection
    SALMON_QUANT(
        SEQKIT_SAMPLE.out.reads,
        SALMON_INDEX.out.index,
        ch_gtf,
        libtype,
    )
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())

    emit:
    reads = SEQKIT_SAMPLE.out.reads // channel: [ val(meta), [ reads ] ]
    index = SALMON_INDEX.out.index // channel: [ index ]
    quant_dir = SALMON_QUANT.out.quant_dir // channel: [ val(meta), [ quant_dir ] ]
    strandness = SALMON_QUANT.out.strandness // channel: [ val(meta), [ strandness ] ]
    quant_log = SALMON_QUANT.out.log // channel: [ val(meta), [ log ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
