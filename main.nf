#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LoaiEletr/rnaseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/LoaiEletr/rnaseq
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RNASEQ } from './workflows/rnaseq'
include { PREPARE_GENOME } from './subworkflows/local/prepare_genome'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { PIPELINE_COMPLETION } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { getGenomeAttribute } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { getContaminantGenome } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta = getGenomeAttribute('genome')
params.transcriptome = getGenomeAttribute('transcript')
params.gtf = getGenomeAttribute('gtf')
params.gtf_isoform = getGenomeAttribute('gtf_isoformswitchanalyzer')
params.bbsplit_index = getGenomeAttribute('bbsplit')
params.sortmerna_index = getGenomeAttribute('sortmerna')
params.hisat2_index = getGenomeAttribute('hisat2')
params.salmon_index = getGenomeAttribute('salmon')
params.kallisto_index = getGenomeAttribute('kallisto')
params.bed = getGenomeAttribute('bed')
params.gff = getGenomeAttribute('gff')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONTAMINANT GENOME FOR BBSplit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Contaminant genome for cross-species detection (e.g., PDX, host-pathogen)
params.contaminant_fasta = getContaminantGenome()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    ch_versions = channel.empty()

    //
    // SUBWORKFLOW : Prepare reference genome files
    //
    PREPARE_GENOME(
        params.fasta,
        params.transcriptome,
        params.gtf,
        params.gtf_isoform,
        params.bed,
        params.gff,
        params.contaminant_fasta,
        params.rrna_db,
        params.rrna_db_type,
        params.aggregation,
        params.bbsplit_index,
        params.sortmerna_index,
        params.hisat2_index,
        params.kallisto_index,
        params.salmon_index,
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden,
    )

    //
    // WORKFLOW: Run main workflow
    //
    LOAIELETR_RNASEQ(
        PIPELINE_INITIALISATION.out.samplesheet,
        ch_versions,
        PREPARE_GENOME.out.transcriptome,
        PREPARE_GENOME.out.fasta_uncompressed,
        PREPARE_GENOME.out.gtf_compressed,
        PREPARE_GENOME.out.gtf_uncompressed,
        PREPARE_GENOME.out.gtf_isoform,
        PREPARE_GENOME.out.rrna_db_fasta,
        PREPARE_GENOME.out.bed,
        PREPARE_GENOME.out.gff,
        PREPARE_GENOME.out.hisat2_index,
        PREPARE_GENOME.out.kallisto_index,
        PREPARE_GENOME.out.salmon_index,
        PREPARE_GENOME.out.sortmerna_index,
        PREPARE_GENOME.out.bbsplit_index,
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        LOAIELETR_RNASEQ.out.multiqc_report,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow LOAIELETR_RNASEQ {
    take:
    ch_samplesheet // channel: [ val(meta), [ reads ] ]
    ch_versions // channel: [ versions.yml ]
    ch_transcriptome // channel: [ transcriptome.fasta ]
    ch_fasta_uncompressed // channel: [ genome.fasta ]
    ch_gtf_compressed // channel: [ genome.gtf.gz ]
    ch_gtf_uncompressed // channel: [ genome.gtf ]
    ch_gtf_isoform // channel: [ isoform.gtf ]
    ch_rrna_db_fasta // channel: [ rrna.fasta ]
    ch_bed // channel: [ genes.bed ]
    ch_gff // channel: [ genes.gff ]
    ch_hisat2_index // channel: [ hisat2_index_files ]
    ch_kallisto_index // channel: [ kallisto_index_files ]
    ch_salmon_index // channel: [ salmon_index_files ]
    ch_sortmerna_index // channel: [ sortmerna_index_files ]
    ch_bbsplit_index // channel: [ bbsplit_index_files ]

    main:

    //
    // WORKFLOW: Run pipeline
    //
    RNASEQ(
        ch_samplesheet,
        ch_versions,
        ch_transcriptome.ifEmpty([]),
        ch_fasta_uncompressed.ifEmpty([]),
        ch_gtf_compressed,
        ch_gtf_uncompressed.ifEmpty([]),
        ch_gtf_isoform.ifEmpty([]),
        ch_rrna_db_fasta.ifEmpty([]),
        ch_bed.ifEmpty([]),
        ch_gff.ifEmpty([]),
        ch_hisat2_index.ifEmpty([]),
        ch_kallisto_index.ifEmpty([]),
        ch_salmon_index.ifEmpty([]),
        ch_sortmerna_index.ifEmpty([]),
        ch_bbsplit_index.ifEmpty([]),
    )

    emit:
    multiqc_report = RNASEQ.out.multiqc_report // channel: /path/to/multiqc_report.html
    versions = RNASEQ.out.versions // channel: [ path(versions.yml) ]
}
