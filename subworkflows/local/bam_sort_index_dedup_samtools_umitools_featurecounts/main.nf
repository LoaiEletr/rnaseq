//
// Process BAM files: sort, index, UMI deduplication, and feature counting
// Performs: SAMtools sort/index, UMI-tools dedup (optional), Subread featureCounts
//

include {
    SAMTOOLS_SORT ;
    SAMTOOLS_SORT as SAMTOOLS_SORT_DEDUP
} from '../../../modules/local/samtools/sort/main.nf'

include {
    SAMTOOLS_INDEX ;
    SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUP
} from '../../../modules/local/samtools/index/main.nf'
include { UMITOOLS_DEDUP } from '../../../modules/local/umitools/dedup/main.nf'
include { SUBREAD_FEATURECOUNTS } from '../../../modules/local/subread/featurecounts/main.nf'

workflow BAM_SORT_INDEX_DEDUP_SAMTOOLS_UMITOOLS_FEATURECOUNTS {
    take:
    ch_bam // channel: [ val(meta), bam ]
    ch_gtf // channel: [ gtf ]
    val_with_umi // boolean: true/false
    val_analysis_method // list: user-selected methods to be executed

    main:

    ch_versions = channel.empty()

    // 1. Sort BAM with SAMtools
    SAMTOOLS_SORT(ch_bam)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())
    ch_umi_bam = SAMTOOLS_SORT.out.sorted_bam

    // 2. Index sorted BAM
    SAMTOOLS_INDEX(ch_umi_bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
    ch_bai = SAMTOOLS_INDEX.out.bai

    // 3. UMI deduplication (optional)
    ch_edit_distance = channel.empty()
    ch_per_umi = channel.empty()
    ch_per_position = channel.empty()
    ch_umi_log = channel.empty()
    if (val_with_umi) {
        UMITOOLS_DEDUP(ch_umi_bam.join(ch_bai))
        ch_versions = ch_versions.mix(UMITOOLS_DEDUP.out.versions.first())
        ch_edit_distance = UMITOOLS_DEDUP.out.tsv_edit_distance
        ch_per_umi = UMITOOLS_DEDUP.out.tsv_per_umi
        ch_per_position = UMITOOLS_DEDUP.out.tsv_umi_per_position
        ch_umi_log = UMITOOLS_DEDUP.out.log

        // 4. Re-sort deduplicated BAM
        SAMTOOLS_SORT_DEDUP(UMITOOLS_DEDUP.out.bam)
        ch_versions = ch_versions.mix(SAMTOOLS_SORT_DEDUP.out.versions.first())
        ch_umi_bam = SAMTOOLS_SORT_DEDUP.out.sorted_bam

        // 5. Re-index deduplicated BAM
        SAMTOOLS_INDEX_DEDUP(ch_umi_bam)
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_DEDUP.out.versions.first())
        ch_bai = SAMTOOLS_INDEX_DEDUP.out.bai
    }

    // 6. Feature counting with Subread
    ch_counts = channel.empty()
    ch_summary = channel.empty()
    if ("DEG" in val_analysis_method) {
        SUBREAD_FEATURECOUNTS(ch_umi_bam, ch_gtf)
        ch_counts = SUBREAD_FEATURECOUNTS.out.counts
        ch_summary = SUBREAD_FEATURECOUNTS.out.summary
        ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())
    }

    emit:
    bam = ch_umi_bam // channel: [ val(meta), [ bam ] ]
    bai = ch_bai // channel: [ val(meta), [ bai ] ]
    edit_distance = ch_edit_distance // channel: [ val(meta), [ tsv_edit_distance ] ]
    per_umi = ch_per_umi // channel: [ val(meta), [ tsv_per_umi ] ]
    per_position = ch_per_position // channel: [ val(meta), [ tsv_umi_per_position ] ]
    umi_log = ch_umi_log // channel: [ val(meta), [ log ] ]
    counts = ch_counts // channel: [ val(meta), [ counts ] ]
    summary = ch_summary // channel: [ val(meta), [ summary ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
