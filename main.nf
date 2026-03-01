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
params.fai = getGenomeAttribute('fai')
params.dict = getGenomeAttribute('dict')
params.interval_list = getGenomeAttribute('interval_list')
params.dbsnp = getGenomeAttribute('dbsnp')
params.dbsnp_tbi = getGenomeAttribute('dbsnp_tbi')
params.knownsites = getGenomeAttribute('knownsites')
params.knownsites_tbi = getGenomeAttribute('knownsites_tbi')
params.snpeff_db = getGenomeAttribute('snpeff_db')
params.snpeff_genome = getGenomeAttribute('snpeff_genome')

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
        params.fai,
        params.dict,
        params.interval_list,
        params.dbsnp,
        params.dbsnp_tbi,
        params.knownsites,
        params.knownsites_tbi,
        params.snpeff_db,
        params.snpeff_genome,
    )

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

    // --- GVC (Germline Variant Calling) Specific Validations ---
    def requested_methods = params.analysis_method ? params.analysis_method.split(",").collect { it.trim() } : []
    if (requested_methods.contains("GVC") && params.aligner == "hisat2") {

        // Handle comma-separated known_sites
        def known_sites_list = params.knownsites ? params.knownsites.split(',').collect { it.trim() }.findAll { it } : []

        // Handle dbSNP Warnings
        if (!params.dbsnp) {
            def manual_species = ["fruitfly", "worm", "monkey"]

            if (manual_species.contains(params.species)) {
                log.warn(
                    """
                ⚠️  MANUAL dbSNP RECOMMENDED: Built-in dbSNP resources are not available for '${params.species}'.
                   When running GVC for this species, it is highly recommended to provide a VCF via: --dbsnp
                """.stripIndent()
                )
            }
            else {
                log.warn(
                    """
                ⚠️  MISSING dbSNP: No dbSNP file provided. HaplotypeCaller will run, but variants
                   will not be labeled with rsIDs. Recommended: --dbsnp 'path/to/dbsnp.vcf'
                """.stripIndent()
                )
            }
        }

        // Ensure Known Sites are provided for BQSR
        if (known_sites_list.isEmpty() && !params.skip_baserecalibration) {
            error(
                """
            ❌ MISSING KNOWN SITES: BQSR requires at least one known variant site file to mask real biology.
               Please either:
                 • Provide sites (comma-separated): --known_sites 'site1.vcf,site2.vcf'
                 • OR skip BQSR using: --skip_baserecalibration
            """.stripIndent()
            )
        }
    }

    //
    // WORKFLOW: Run main workflow
    //
    LOAIELETR_RNASEQ(
        PIPELINE_INITIALISATION.out.samplesheet,
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
        PREPARE_GENOME.out.fai,
        PREPARE_GENOME.out.dict,
        PREPARE_GENOME.out.intervallist,
        PREPARE_GENOME.out.intervals_split,
        PREPARE_GENOME.out.knownsites,
        PREPARE_GENOME.out.knownsites_tbi,
        PREPARE_GENOME.out.dbsnp,
        PREPARE_GENOME.out.dbsnp_tbi,
        PREPARE_GENOME.out.snpeff_db,
        PREPARE_GENOME.out.snpeff_genome,
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
    ch_fai // channel: [ val(meta), [ fai ] ]
    ch_dict // channel: [ val(meta), [ dict ] ]
    ch_intervallist // channel: [ val(meta), [ interval_list ] ]
    ch_intervals_split // channel: [ val(meta), [ interval_list ] ]
    ch_knownsites // channel: [ val(meta), [ knownsites ] ]
    ch_knownsites_tbi // channel: [ val(meta), [ knownsites_tbi ] ]
    ch_dbsnp // channel: [ val(meta), [ dbsnp ] ]
    ch_dbsnp_tbi // channel: [ val(meta), [ dbsnp_tbi ] ]
    ch_snpeff_db // channel: [ val(meta), [ snpeff_db ] ]
    ch_snpeff_genome // channel: [ val(meta), [ snpeff_genome ] ]

    main:

    //
    // WORKFLOW: Run pipeline
    //
    RNASEQ(
        ch_samplesheet,
        ch_transcriptome.ifEmpty([[:], []]),
        ch_fasta_uncompressed.ifEmpty([[:], []]),
        ch_gtf_compressed,
        ch_gtf_uncompressed.ifEmpty([[:], []]),
        ch_gtf_isoform.ifEmpty([[:], []]),
        ch_rrna_db_fasta.ifEmpty([[:], []]),
        ch_bed.ifEmpty([[:], []]),
        ch_gff.ifEmpty([[:], []]),
        ch_hisat2_index.ifEmpty([[:], []]),
        ch_kallisto_index.ifEmpty([[:], []]),
        ch_salmon_index.ifEmpty([[:], []]),
        ch_sortmerna_index.ifEmpty([[:], []]),
        ch_bbsplit_index.ifEmpty([[:], []]),
        ch_fai.ifEmpty([[:], []]),
        ch_dict.ifEmpty([[:], []]),
        ch_intervallist.ifEmpty([[:], []]),
        ch_intervals_split.ifEmpty([[:], []]),
        ch_knownsites.ifEmpty([[:], []]),
        ch_knownsites_tbi.ifEmpty([[:], []]),
        ch_dbsnp.ifEmpty([[:], []]),
        ch_dbsnp_tbi.ifEmpty([[:], []]),
        ch_snpeff_db.ifEmpty([[:], []]),
        ch_snpeff_genome.ifEmpty([[:], []]),
    )

    emit:
    multiqc_report = RNASEQ.out.multiqc_report // channel: /path/to/multiqc_report.html
}
