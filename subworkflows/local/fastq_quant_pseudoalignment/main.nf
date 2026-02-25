//
// Quantify with pseudoaligner (Kallisto or Salmon)
// Performs: Transcriptome quantification with user-selected pseudoaligner
//

include { KALLISTO_QUANT } from '../../../modules/local/kallisto/quant/main.nf'
include { SALMON_QUANT } from '../../../modules/local/salmon/quant/main.nf'

workflow FASTQ_QUANT_PSEUDOALIGNMENT {
    take:
    ch_reads // channel: [ val(meta), reads ]
    ch_kallisto_index // channel: [ val(meta), kallisto_index ]
    ch_salmon_index // channel: [ val(meta), salmon_index ]
    ch_gtf // channel: [ val(meta), gtf ]
    val_fragment_length // integer: Fragment length for Kallisto
    val_fragment_length_sd // integer: Fragment length standard deviation for Kallisto
    val_pseudoaligner_method // string: 'kallisto' or 'salmon'

    main:

    ch_quant_dir = channel.empty()
    ch_quant_log = channel.empty()

    if (val_pseudoaligner_method == "kallisto") {

        KALLISTO_QUANT(
            ch_reads,
            ch_kallisto_index,
            ch_gtf,
            val_fragment_length,
            val_fragment_length_sd,
        )
        ch_quant_dir = KALLISTO_QUANT.out.quant_dir
        ch_quant_log = KALLISTO_QUANT.out.log
    }
    else if (val_pseudoaligner_method == "salmon") {

        SALMON_QUANT(
            ch_reads,
            ch_salmon_index,
            ch_gtf,
        )
        ch_quant_dir = SALMON_QUANT.out.quant_dir
    }
    else {
        error("Invalid pseudoaligner_method: ${val_pseudoaligner_method}. Must be one of: 'salmon' or 'kallisto'")
    }

    emit:
    quant_dir = ch_quant_dir // channel: [ val(meta), [ quant_dir ] ]
    quant_log = ch_quant_log // channel: [ val(meta), [ quant_log ] ] (Kallisto only)
}
