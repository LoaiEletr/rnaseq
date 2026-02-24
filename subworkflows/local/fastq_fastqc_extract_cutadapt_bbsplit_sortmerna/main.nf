//
// FASTQ quality control, UMI processing, adapter trimming, contamination removal, and rRNA filtering
// Performs: FASTQC (raw/trimmed), UMI extraction, Cutadapt trimming, BBMap bbsplit contamination removal, SortMeRNA rRNA filtering
//

include {
    FASTQC as FASTQC_RAW ;
    FASTQC as FASTQC_TRIM
} from '../../../modules/local/fastqc/main.nf'
include { CUTADAPT } from '../../../modules/local/cutadapt/main.nf'
include { UMITOOLS_EXTRACT } from '../../../modules/local/umitools/extract/main.nf'
include { BBMAP_BBSPLIT } from '../../../modules/local/bbmap/bbsplit/main.nf'
include { SORTMERNA } from '../../../modules/local/sortmerna/main.nf'

workflow FASTQ_FASTQC_EXTRACT_CUTADAPT_BBSPLIT_SORTMERNA {
    take:
    ch_reads // channel: [ val(meta), reads ]
    ch_adapter // channel: [ val(meta), adapter ]
    ch_rrna_db // channel: [ val(meta), rrna_db ]
    ch_bbsplit_index // channel: [ val(meta), bbsplit_index ]
    ch_sortmerna_index // channel: [ val(meta), sortmerna_index ]
    val_rrna_db_type // string: default, fast, or sensitive
    val_with_umi // boolean: true/false
    val_skip_umi_extract // boolean: true/false
    val_skip_fastqc // boolean: true/false
    val_skip_cutadapt // boolean: true/false
    val_skip_bbsplit // boolean: true/false
    val_skip_sortmerna // boolean: true/false

    main:

    ch_fastqc_raw_html = channel.empty()
    ch_fastqc_raw_zip = channel.empty()

    // 1. Raw reads quality control (optional)
    if (!val_skip_fastqc) {
        FASTQC_RAW(ch_reads)
        ch_fastqc_raw_html = FASTQC_RAW.out.html
        ch_fastqc_raw_zip = FASTQC_RAW.out.zip
    }

    ch_umi_reads = ch_reads
    ch_umi_log = channel.empty()

    // 2. UMI extraction (optional, requires with_umi=true)
    if (val_with_umi && !val_skip_umi_extract) {
        UMITOOLS_EXTRACT(ch_umi_reads)
        ch_umi_reads = UMITOOLS_EXTRACT.out.reads
        ch_umi_log = UMITOOLS_EXTRACT.out.log
    }

    ch_trim_reads = ch_umi_reads
    ch_trim_json = channel.empty()

    // 3. Adapter trimming (optional)
    if (!val_skip_cutadapt) {
        CUTADAPT(ch_trim_reads, ch_adapter)
        ch_trim_reads = CUTADAPT.out.reads
        ch_trim_json = CUTADAPT.out.json
    }

    ch_fastqc_trim_html = channel.empty()
    ch_fastqc_trim_zip = channel.empty()

    // 4. Trimmed reads quality control (optional, requires trimming)
    if (!val_skip_fastqc && !val_skip_cutadapt) {
        FASTQC_TRIM(ch_trim_reads)
        ch_fastqc_trim_html = FASTQC_TRIM.out.html
        ch_fastqc_trim_zip = FASTQC_TRIM.out.zip
    }

    ch_bbsplit_reads = ch_trim_reads
    ch_bbsplit_stats = channel.empty()

    // 5. Contamination removal with BBMap bbsplit (optional)
    if (!val_skip_bbsplit) {
        BBMAP_BBSPLIT(
            ch_bbsplit_reads,
            ch_bbsplit_index,
            [[:], []],
            [[:], []],
            false,
        )
        ch_bbsplit_reads = BBMAP_BBSPLIT.out.reads
        ch_bbsplit_stats = BBMAP_BBSPLIT.out.stats
    }

    ch_nonrrna_reads = ch_bbsplit_reads
    ch_rrna_log = channel.empty()

    // 6. rRNA filtering with SortMeRNA (optional)
    if (!val_skip_sortmerna) {
        // Filter rRNA reads
        SORTMERNA(
            ch_nonrrna_reads,
            ch_sortmerna_index,
            ch_rrna_db,
            val_rrna_db_type,
            false,
        )
        ch_nonrrna_reads = SORTMERNA.out.reads
        ch_rrna_log = SORTMERNA.out.log
    }

    emit:
    fastqc_raw_html = ch_fastqc_raw_html // channel: [ val(meta), [ html ] ]
    fastqc_raw_zip = ch_fastqc_raw_zip // channel: [ val(meta), [ zip ] ]
    umi_log = ch_umi_log // channel: [ val(meta), [ log ] ]
    trim_json = ch_trim_json // channel: [ val(meta), [ json ] ]
    fastqc_trim_html = ch_fastqc_trim_html // channel: [ val(meta), [ html ] ]
    fastqc_trim_zip = ch_fastqc_trim_zip // channel: [ val(meta), [ zip ] ]
    bbsplit_stats = ch_bbsplit_stats // channel: [ val(meta), [ stats ] ]
    reads = ch_nonrrna_reads // channel: [ val(meta), [ reads ] ]
    rrna_log = ch_rrna_log // channel: [ val(meta), [ log ] ]
}
