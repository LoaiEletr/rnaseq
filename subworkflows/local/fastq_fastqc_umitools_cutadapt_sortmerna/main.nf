//
// FASTQ quality control, UMI extraction, adapter trimming, and rRNA filtering workflow
// Performs: FASTQC (raw/trimmed), UMI extraction, Cutadapt trimming, SortMeRNA rRNA filtering
//

include { FASTQC as FASTQC_RAW } from '../../../modules/local/fastqc/main.nf'
include { FASTQC as FASTQC_TRIM } from '../../../modules/local/fastqc/main.nf'
include { CUTADAPT } from '../../../modules/local/cutadapt/main.nf'
include { UMITOOLS_EXTRACT } from '../../../modules/local/umitools/extract/main.nf'
include {
    SORTMERNA as SORTMERNA_INDEX ;
    SORTMERNA as SORTMERNA_ALIGN
} from '../../../modules/local/sortmerna/main.nf'

workflow FASTQ_FASTQC_UMITOOLS_CUTADAPT_SORTMERNA {
    take:
    ch_reads // channel: [ val(meta), reads ]
    ch_adapter // channel: [ adapter ]
    ch_rrna_db // channel: [ rrna_db ]
    val_rrna_db_type // string: default, fast, or sensitive
    val_lib_kit // string: quantseq, corall, or takara
    val_with_umi // boolean: true/false
    val_skip_umi_extract // boolean: true/false
    val_skip_fastqc // boolean: true/false
    val_skip_cutadapt // boolean: true/false
    val_skip_sortmerna // boolean: true/false

    main:

    ch_versions = channel.empty()
    ch_fastqc_raw_html = channel.empty()
    ch_fastqc_raw_zip = channel.empty()

    // 1. Raw reads quality control (optional)
    if (!val_skip_fastqc) {
        FASTQC_RAW(ch_reads)
        ch_fastqc_raw_html = FASTQC_RAW.out.html
        ch_fastqc_raw_zip = FASTQC_RAW.out.zip
        ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())
    }

    ch_umi_reads = ch_reads
    ch_umi_log = channel.empty()

    // 2. UMI extraction (optional, requires with_umi=true)
    if (val_with_umi && !val_skip_umi_extract) {
        UMITOOLS_EXTRACT(ch_umi_reads, val_lib_kit)
        ch_umi_reads = UMITOOLS_EXTRACT.out.reads
        ch_umi_log = UMITOOLS_EXTRACT.out.log
        ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions.first())
    }

    ch_trim_reads = ch_umi_reads
    ch_trim_json = channel.empty()

    // 3. Adapter trimming (optional)
    if (!val_skip_cutadapt) {
        CUTADAPT(ch_trim_reads, ch_adapter)
        ch_trim_reads = CUTADAPT.out.reads
        ch_trim_json = CUTADAPT.out.json
        ch_versions = ch_versions.mix(CUTADAPT.out.versions.first())
    }

    ch_fastqc_trim_html = channel.empty()
    ch_fastqc_trim_zip = channel.empty()

    // 4. Trimmed reads quality control (optional, requires trimming)
    if (!val_skip_fastqc && !val_skip_cutadapt) {
        FASTQC_TRIM(ch_trim_reads)
        ch_fastqc_trim_html = FASTQC_TRIM.out.html
        ch_fastqc_trim_zip = FASTQC_TRIM.out.zip
        ch_versions = ch_versions.mix(FASTQC_TRIM.out.versions.first())
    }

    ch_nonrrna_reads = ch_trim_reads
    ch_nonrrna_log = channel.empty()
    ch_index = channel.empty()

    // 5. rRNA filtering with SortMeRNA (optional)
    if (!val_skip_sortmerna) {
        // Build rRNA index
        SORTMERNA_INDEX(
            ch_reads,
            [],
            ch_rrna_db,
            val_rrna_db_type,
            true,
        )
        ch_index = SORTMERNA_INDEX.out.index
        ch_versions = ch_versions.mix(SORTMERNA_INDEX.out.versions)

        // Filter rRNA reads
        SORTMERNA_ALIGN(
            ch_reads,
            ch_index,
            ch_rrna_db,
            val_rrna_db_type,
            false,
        )
        ch_nonrrna_reads = SORTMERNA_ALIGN.out.reads
        ch_nonrrna_log = SORTMERNA_ALIGN.out.log
        ch_versions = ch_versions.mix(SORTMERNA_ALIGN.out.versions.first())
    }

    emit:
    fastqc_raw_html = ch_fastqc_raw_html // channel: [ val(meta), [ html ] ]
    fastqc_raw_zip = ch_fastqc_raw_zip // channel: [ val(meta), [ zip ] ]
    umi_log = ch_umi_log // channel: [ val(meta), [ log ] ]
    trim_json = ch_trim_json // channel: [ val(meta), [ json ] ]
    fastqc_trim_html = ch_fastqc_trim_html // channel: [ val(meta), [ html ] ]
    fastqc_trim_zip = ch_fastqc_trim_zip // channel: [ val(meta), [ zip ] ]
    index = ch_index // channel: [ index ]
    reads = ch_nonrrna_reads // channel: [ val(meta), [ reads ] ]
    nonrrna_log = ch_nonrrna_log // channel: [ val(meta), [ log ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
