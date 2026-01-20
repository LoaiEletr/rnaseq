//
// Index transcriptome and quantify with pseudoaligner (Kallisto or Salmon)
// Performs: Transcriptome indexing and quantification with user-selected pseudoaligner
//

include { KALLISTO_INDEX } from '../../../modules/local/kallisto/index/main.nf'
include { SALMON_INDEX } from '../../../modules/local/salmon/index/main.nf'
include { KALLISTO_QUANT } from '../../../modules/local/kallisto/quant/main.nf'
include { SALMON_QUANT } from '../../../modules/local/salmon/quant/main.nf'

workflow FASTQ_INDEX_QUANT_PSEUDOALIGNMENT {
    take:
    reads // channel: [ val(meta), reads ]
    ch_transcriptome // channel: [ transcriptome ]
    ch_gtf // channel: [ gtf ]
    val_bootstrap_count // integer: Bootstrap count for Kallisto
    val_fragment_length // integer: Fragment length for Kallisto
    val_fragment_length_sd // integer: Fragment length standard deviation for Kallisto
    val_pseudoaligner_method // string: 'kallisto' or 'salmon'

    main:

    ch_versions = channel.empty()
    ch_index = channel.empty()
    ch_quant_dir = channel.empty()
    ch_json_info = channel.empty()
    ch_quant_log = channel.empty()

    if (val_pseudoaligner_method == "kallisto") {
        KALLISTO_INDEX(ch_transcriptome)
        ch_versions = ch_versions.mix(KALLISTO_INDEX.out.versions)
        ch_index = KALLISTO_INDEX.out.index

        KALLISTO_QUANT(
            reads,
            ch_index,
            ch_gtf,
            val_bootstrap_count,
            val_fragment_length,
            val_fragment_length_sd,
        )
        ch_versions = ch_versions.mix(KALLISTO_QUANT.out.versions.first())
        ch_quant_dir = KALLISTO_QUANT.out.quant_dir
        ch_json_info = KALLISTO_QUANT.out.json_info
        ch_quant_log = KALLISTO_QUANT.out.log
    }
    else if (val_pseudoaligner_method == "salmon") {
        SALMON_INDEX(ch_transcriptome)
        ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)
        ch_index = SALMON_INDEX.out.index

        SALMON_QUANT(
            reads,
            ch_index,
            ch_gtf,
            [],
        )
        ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())
        ch_quant_dir = SALMON_QUANT.out.quant_dir
        ch_quant_log = SALMON_QUANT.out.log
    }
    else {
        error("Invalid pseudoaligner_method: ${val_pseudoaligner_method}. Must be one of: 'salmon' or 'kallisto'")
    }

    emit:
    index = ch_index // channel: [ index ]
    quant_dir = ch_quant_dir // channel: [ val(meta), [ quant_dir ] ]
    json_info = ch_json_info // channel: [ val(meta), [ json_info ] ] (Kallisto only)
    quant_log = ch_quant_log // channel: [ val(meta), [ quant_log ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
