//
// Read deduplication, RNA-seq transcriptomic splitting, and Base Quality Score Recalibration
//
// Performs: Picard MarkDuplicates, GATK4 SplitNCigarReads, GATK4 BaseRecalibrator, and GATK4 ApplyBQSR
//

include { PICARD_MARKDUPLICATES } from '../../../modules/local/picard/markduplicates'
include { GATK4_SPLITNCIGARREADS } from '../../../modules/local/gatk4/splitncigarreads'
include { SAMTOOLS_MERGE } from '../../../modules/local/samtools/merge'
include {
    SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUP ;
    SAMTOOLS_INDEX as SAMTOOLS_INDEX_MERGED ;
    SAMTOOLS_INDEX as SAMTOOLS_INDEX_BQSR
} from '../../../modules/local/samtools/index'
include { GATK4_BASERECALIBRATOR } from '../../../modules/local/gatk4/baserecalibrator'
include { GATK4_APPLYBQSR } from '../../../modules/local/gatk4/applybqsr'

workflow BAM_DEDUP_RECALIBRATE_PICARD_GATK {
    take:
    ch_bam // channel: [ val(meta), [ bam ] ]
    ch_bai // channel: [ val(meta), [ bai ] ]
    ch_fasta // channel: [ val(meta), [ fasta ] ]
    ch_fai // channel: [ val(meta), [ fai ] ]
    ch_dict // channel: [ val(meta), [ dict ] ]
    ch_intervals // channel: [ val(meta), [ intervals ] ]
    ch_known_sites // channel: [ val(meta), [ vcf ] ]
    ch_known_sites_tbi // channel: [ val(meta), [ tbi ] ]
    skip_picard_markduplicates // boolean: true/false
    skip_baserecalibration // boolean: true/false

    main:

    // 1. Mark duplicate reads using Picard
    ch_picard_bam = ch_bam
    ch_picard_bai = ch_bai
    ch_picard_metrics = channel.empty()

    if (!skip_picard_markduplicates) {
        PICARD_MARKDUPLICATES(ch_bam)
        ch_picard_metrics = PICARD_MARKDUPLICATES.out.metrics
        ch_picard_bam = PICARD_MARKDUPLICATES.out.bam

        SAMTOOLS_INDEX_DEDUP(ch_picard_bam)
        ch_picard_bai = SAMTOOLS_INDEX_DEDUP.out.bai
    }

    // 2. Split reads with Ns in CIGAR (RNA-seq specific processing)
    GATK4_SPLITNCIGARREADS(
        ch_picard_bam.join(ch_picard_bai).combine(ch_intervals.map { it[1] }),
        ch_fasta,
        ch_fai,
        ch_dict,
    )

    // 3. Merge and index split BAM files
    SAMTOOLS_MERGE(GATK4_SPLITNCIGARREADS.out.bam.groupTuple())
    ch_merged_bam = SAMTOOLS_MERGE.out.bam

    SAMTOOLS_INDEX_MERGED(ch_merged_bam)
    ch_merged_bai = SAMTOOLS_INDEX_MERGED.out.bai

    // 4. Generate BQSR recalibration table
    ch_recalibrate_bam = ch_merged_bam
    ch_recalibrate_bai = ch_merged_bai

    if (!skip_baserecalibration) {
        GATK4_BASERECALIBRATOR(
            ch_recalibrate_bam.join(ch_recalibrate_bai).combine(ch_intervals.map { it[1] }),
            ch_fasta,
            ch_fai,
            ch_dict,
            ch_known_sites,
            ch_known_sites_tbi,
        )

        // 5. Apply base quality score recalibration
        GATK4_APPLYBQSR(
            ch_recalibrate_bam.join(ch_recalibrate_bai).join(GATK4_BASERECALIBRATOR.out.table).combine(ch_intervals.map { it[1] }),
            ch_fasta,
            ch_fai,
            ch_dict,
        )

        ch_recalibrate_bam = GATK4_APPLYBQSR.out.bam
        SAMTOOLS_INDEX_BQSR(ch_recalibrate_bam)
        ch_recalibrate_bai = SAMTOOLS_INDEX_BQSR.out.bai
    }

    emit:
    bam = ch_recalibrate_bam // channel: [ val(meta), [ bam ] ]
    bai = ch_recalibrate_bai // channel: [ val(meta), [ bai ] ]
    picard_metrics = ch_picard_metrics // channel: [ val(meta), [ metrics ] ]
}
