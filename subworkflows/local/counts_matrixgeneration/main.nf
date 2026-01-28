//
// Generate gene count matrix from quantification results
// Performs: Transcript-to-gene mapping + tximport (for Salmon/Kallisto) or count merging (for featureCounts)
//

include { TX2GENE } from '../../../modules/local/tx2gene/main.nf'
include { TXIMPORT } from '../../../modules/local/tximport/main.nf'
include { MERGE_COUNTS } from '../../../modules/local/merge_counts/main.nf'

workflow COUNTS_MATRIXGENERATION {
    take:
    ch_counts // channel: [ quant_results ]
    val_quant_method // string: 'salmon', 'kallisto', or 'featurecounts'
    val_species_name // string: Species name for transcript-to-gene mapping

    main:

    ch_versions = channel.empty()
    ch_tx2gene_csv = channel.empty()
    ch_gene_table_rds = channel.empty()

    // 1. For Salmon/Kallisto: Transcript-to-gene mapping + tximport
    if (val_quant_method == "salmon" || val_quant_method == "kallisto") {
        TX2GENE(val_species_name)
        ch_tx2gene_csv = TX2GENE.out.tsv
        ch_versions = ch_versions.mix(TX2GENE.out.versions)

        TXIMPORT(
            ch_counts,
            val_quant_method,
            TX2GENE.out.tsv,
        )
        ch_gene_table_rds = TXIMPORT.out.rds
        ch_versions = ch_versions.mix(TXIMPORT.out.versions)
    }
    else if (val_quant_method == "featurecounts") {
        MERGE_COUNTS(
            ch_counts
        )
        ch_gene_table_rds = MERGE_COUNTS.out.counts_table
        ch_versions = ch_versions.mix(MERGE_COUNTS.out.versions)
    }
    else {
        error("Invalid quant_method: ${val_quant_method}. Must be one of: 'salmon', 'kallisto', or 'featurecounts'")
    }

    emit:
    gene_table_rds = ch_gene_table_rds // channel: [ rds ]
    versions = ch_versions // channel: [ versions.yml ]
}
