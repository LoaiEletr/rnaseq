//
// Differential splicing analysis with rMATS
// Performs: BAM list creation, read length calculation, rMATS preparation and post-processing,
//           significant event filtering, and Sashimi plot generation
//

include { MAKE_BAMLIST } from '../../../modules/local/make_bamlist/main.nf'
include { SEQKIT_STATS } from '../../../modules/local/seqkit/stats/main.nf'
include {
    RMATS as RMATS_PREP ;
    RMATS as RMATS_POST
} from '../../../modules/local/rmats/main.nf'
include { FILTER_SIGEVENTS } from '../../../modules/local/filter_sigevents/main.nf'
include { RMATS2SASHIMIPLOT } from '../../../modules/local/rmats2sashimiplot/main.nf'

workflow BAM_BAMLIST_AVGLEN_RMATS {
    take:
    ch_bam // channel: [ val(meta), bams ]
    ch_fastq_subsample // channel: [ val(meta), fastq ]
    ch_gtf // channel: [ val(meta), gtf ]
    val_statoff // boolean: true/false: Turn off statistical analysis
    val_novelss // boolean: true/false: Detect novel splicing sites
    val_allowclipping // boolean: true/false: Allow clipped alignments
    val_individualcounts // boolean: true/false: Output individual read counts
    val_event_types // string: Event types to analyze (e.g., 'SE,RI')
    val_fdr_cutoff // float: FDR cutoff for filtering
    val_delta_psi // float: ΔPSI cutoff for filtering

    main:

    ch_versions = channel.empty()

    // 1. Create BAM list for rMATS input
    MAKE_BAMLIST(
        ch_bam.map { meta, bams ->
            [[condition: meta.condition, lib_type: meta.lib_type], bams]
        }.groupTuple()
    )
    ch_versions = ch_versions.mix(MAKE_BAMLIST.out.versions.first())

    // 2. Calculate average read length from subsampled FASTQ
    SEQKIT_STATS(ch_fastq_subsample)
    ch_versions = ch_versions.mix(SEQKIT_STATS.out.versions.first())

    // 3. rMATS preparation step
    RMATS_PREP(
        ch_bam.map { meta, bams ->
            [[single_end: meta.single_end, lib_type: meta.lib_type], bams]
        }.groupTuple(),
        MAKE_BAMLIST.out.txt.map { meta, bams -> bams }.collect(),
        ch_gtf,
        'prep',
        SEQKIT_STATS.out.avg_length,
        val_statoff,
        val_novelss,
        val_allowclipping,
        val_individualcounts,
        [],
    )
    ch_versions = ch_versions.mix(RMATS_PREP.out.versions)

    // 4. rMATS post-processing step
    RMATS_POST(
        ch_bam.map { meta, bams ->
            [[single_end: meta.single_end, lib_type: meta.lib_type], bams]
        }.groupTuple(),
        MAKE_BAMLIST.out.txt.map { meta, bams -> bams }.collect(),
        ch_gtf,
        'post',
        SEQKIT_STATS.out.avg_length,
        val_statoff,
        val_novelss,
        val_allowclipping,
        val_individualcounts,
        RMATS_PREP.out.prep_tmp,
    )
    ch_versions = ch_versions.mix(RMATS_POST.out.versions)

    // 5. Filter significant splicing events
    FILTER_SIGEVENTS(
        RMATS_POST.out.post_output,
        val_event_types,
        val_fdr_cutoff,
        val_delta_psi,
    )
    ch_versions = ch_versions.mix(FILTER_SIGEVENTS.out.versions.first())

    // 6. Generate Sashimi plots for significant events
    RMATS2SASHIMIPLOT(
        FILTER_SIGEVENTS.out.sig_rmats.combine(
            MAKE_BAMLIST.out.txt.map { meta, bams -> bams }.collect()
        ),
        ch_bam.map { meta, bams -> bams }.collect(),
    )
    ch_versions = ch_versions.mix(RMATS2SASHIMIPLOT.out.versions.first())

    emit:
    rmats_prep_tmp = RMATS_PREP.out.prep_tmp // channel: [ prep_tmp ]
    rmats_output = RMATS_POST.out.post_output // channel: [ post_output ]
    sigevents_txt = FILTER_SIGEVENTS.out.sig_rmats // channel: [ val(event_type), sig_rmats ]
    rmats2sashimi_output = RMATS2SASHIMIPLOT.out.sashimi_output // channel: [ sashimi_output ]
    versions = ch_versions // channel: [ versions.yml ]
}
