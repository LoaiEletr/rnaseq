//
// Comprehensive BAM quality control with RSeQC modules
// Performs multiple RSeQC analyses based on selected modules
//

include { RSEQC_BAMSTAT } from '../../../modules/local/rseqc/bamstat/main.nf'
include { RSEQC_GENEBODYCOVERAGE } from '../../../modules/local/rseqc/genebodycoverage/main.nf'
include { RSEQC_INFEREXPERIMENT } from '../../../modules/local/rseqc/inferexperiment/main.nf'
include { RSEQC_INNERDISTANCE } from '../../../modules/local/rseqc/innerdistance/main.nf'
include { RSEQC_JUNCTIONANNOTATION } from '../../../modules/local/rseqc/junctionannotation/main.nf'
include { RSEQC_READDISTRIBUTION } from '../../../modules/local/rseqc/readdistribution/main.nf'
include { RSEQC_READDUPLICATION } from '../../../modules/local/rseqc/readduplication/main.nf'
include { RSEQC_TIN } from '../../../modules/local/rseqc/tin/main.nf'

workflow BAM_RSEQC {
    take:
    rseqc_modules // channel: [ modules ] or list of module names
    ch_bam // channel: [ val(meta), bam ]
    ch_bed // channel: [ bed ]
    ch_bai // channel: [ val(meta), bai ]

    main:

    ch_versions = channel.empty()
    ch_bamstat = channel.empty()

    if ("bam_stat" in rseqc_modules) {
        RSEQC_BAMSTAT(ch_bam)
        ch_bamstat = RSEQC_BAMSTAT.out.txt
        ch_versions = ch_versions.mix(RSEQC_BAMSTAT.out.versions.first())
    }

    ch_genebodycoverage_r = channel.empty()
    ch_genebodycoverage_pdf = channel.empty()
    ch_genebodycoverage_txt = channel.empty()
    ch_genebodycoverage_log = channel.empty()

    if ("genebody_coverage" in rseqc_modules) {
        RSEQC_GENEBODYCOVERAGE(
            ch_bam.map { meta, bams -> bams }.collect(),
            ch_bai.map { meta, bais -> bais }.collect(),
            ch_bed,
        )
        ch_genebodycoverage_r = RSEQC_GENEBODYCOVERAGE.out.rscript
        ch_genebodycoverage_pdf = RSEQC_GENEBODYCOVERAGE.out.pdf
        ch_genebodycoverage_txt = RSEQC_GENEBODYCOVERAGE.out.txt
        ch_genebodycoverage_log = RSEQC_GENEBODYCOVERAGE.out.log
        ch_versions = ch_versions.mix(RSEQC_GENEBODYCOVERAGE.out.versions)
    }

    ch_inferexperiment = channel.empty()

    if ("infer_experiment" in rseqc_modules) {
        RSEQC_INFEREXPERIMENT(ch_bam, ch_bed)
        ch_inferexperiment = RSEQC_INFEREXPERIMENT.out.txt
        ch_versions = ch_versions.mix(RSEQC_INFEREXPERIMENT.out.versions.first())
    }

    ch_innerdistance_distance = channel.empty()
    ch_innerdistance_freq = channel.empty()
    ch_innerdistance_pdf = channel.empty()
    ch_innerdistance_r = channel.empty()

    if ("inner_distance" in rseqc_modules) {
        RSEQC_INNERDISTANCE(ch_bam, ch_bed)
        ch_innerdistance_distance = RSEQC_INNERDISTANCE.out.distance
        ch_innerdistance_freq = RSEQC_INNERDISTANCE.out.freq
        ch_innerdistance_pdf = RSEQC_INNERDISTANCE.out.pdf
        ch_innerdistance_r = RSEQC_INNERDISTANCE.out.rscript
        ch_versions = ch_versions.mix(RSEQC_INNERDISTANCE.out.versions.first())
    }

    ch_junctionannotation_xls = channel.empty()
    ch_junctionannotation_r = channel.empty()
    ch_junctionannotation_bed = channel.empty()
    ch_junctionannotation_interact = channel.empty()
    ch_junctionannotation_pdf = channel.empty()
    ch_junctionannotation_events_pdf = channel.empty()

    if ("junction_annotation" in rseqc_modules) {
        RSEQC_JUNCTIONANNOTATION(ch_bam, ch_bed)
        ch_junctionannotation_xls = RSEQC_JUNCTIONANNOTATION.out.xls
        ch_junctionannotation_r = RSEQC_JUNCTIONANNOTATION.out.rscript
        ch_junctionannotation_bed = RSEQC_JUNCTIONANNOTATION.out.junction_bed
        ch_junctionannotation_interact = RSEQC_JUNCTIONANNOTATION.out.interact_bed
        ch_junctionannotation_pdf = RSEQC_JUNCTIONANNOTATION.out.junction_pdf
        ch_junctionannotation_events_pdf = RSEQC_JUNCTIONANNOTATION.out.events_pdf
        ch_versions = ch_versions.mix(RSEQC_JUNCTIONANNOTATION.out.versions.first())
    }

    ch_readdistribution = channel.empty()

    if ("read_distribution" in rseqc_modules) {
        RSEQC_READDISTRIBUTION(ch_bam, ch_bed)
        ch_readdistribution = RSEQC_READDISTRIBUTION.out.txt
        // Assumed output name
        ch_versions = ch_versions.mix(RSEQC_READDISTRIBUTION.out.versions.first())
    }

    ch_readduplication_pdf = channel.empty()
    ch_readduplication_r = channel.empty()
    ch_readduplication_seq = channel.empty()
    ch_readduplication_pos = channel.empty()

    if ("read_duplication" in rseqc_modules) {
        RSEQC_READDUPLICATION(ch_bam)
        ch_readduplication_pdf = RSEQC_READDUPLICATION.out.pdf
        ch_readduplication_r = RSEQC_READDUPLICATION.out.rscript
        ch_readduplication_seq = RSEQC_READDUPLICATION.out.seq_duplicaterate
        ch_readduplication_pos = RSEQC_READDUPLICATION.out.pos_duplicaterate
        ch_versions = ch_versions.mix(RSEQC_READDUPLICATION.out.versions.first())
    }

    ch_tin_txt = channel.empty()
    ch_tin_xls = channel.empty()

    if ("tin" in rseqc_modules) {
        RSEQC_TIN(ch_bam.join(ch_bai), ch_bed)
        ch_tin_txt = RSEQC_TIN.out.txt
        ch_tin_xls = RSEQC_TIN.out.xls
        ch_versions = ch_versions.mix(RSEQC_TIN.out.versions.first())
    }

    emit:
    bamstat = ch_bamstat // channel: [ val(meta), txt ]
    genebodycoverage_r = ch_genebodycoverage_r // channel: [ rscript ]
    genebodycoverage_pdf = ch_genebodycoverage_pdf // channel: [ pdf ]
    genebodycoverage_txt = ch_genebodycoverage_txt // channel: [ txt ]
    genebodycoverage_log = ch_genebodycoverage_log // channel: [ log ]
    inferexperiment = ch_inferexperiment // channel: [ val(meta), [ txt ] ]
    innerdistance_distance = ch_innerdistance_distance // channel: [ val(meta), [ distance ] ]
    innerdistance_freq = ch_innerdistance_freq // channel: [ val(meta), [ freq ] ]
    innerdistance_pdf = ch_innerdistance_pdf // channel: [ val(meta), [ pdf ] ]
    innerdistance_r = ch_innerdistance_r // channel: [ val(meta), [ rscript ] ]
    junctionannotation_xls = ch_junctionannotation_xls // channel: [ val(meta), [ xls ] ]
    junctionannotation_r = ch_junctionannotation_r // channel: [ val(meta), [ rscript ] ]
    junctionannotation_bed = ch_junctionannotation_bed // channel: [ val(meta), [ junction_bed ] ]
    junctionannotation_interact = ch_junctionannotation_interact // channel: [ val(meta), [ interact_bed ] ]
    junctionannotation_pdf = ch_junctionannotation_pdf // channel: [ val(meta), [ junction_pdf ] ]
    junctionannotation_events_pdf = ch_junctionannotation_events_pdf // channel: [ val(meta), [ events_pdf ] ]
    readdistribution = ch_readdistribution // channel: [ val(meta), [ txt ] ]
    readduplication_pdf = ch_readduplication_pdf // channel: [ val(meta), [ pdf ] ]
    readduplication_r = ch_readduplication_r // channel: [ val(meta), [ rscript ] ]
    readduplication_seq = ch_readduplication_seq // channel: [ val(meta), [ seq_duplicaterate ] ]
    readduplication_pos = ch_readduplication_pos // channel: [ val(meta), [ pos_duplicaterate ] ]
    tin_txt = ch_tin_txt // channel: [ val(meta), [ txt ] ]
    tin_xls = ch_tin_xls // channel: [ val(meta), [ xls ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
