/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES: Local hisat2, seqkit, and reporting
//
include { HISAT2_ALIGN } from '../modules/local/hisat2/align'
include { SEQKIT_SAMPLE } from '../modules/local/seqkit/sample'
include { BCFTOOLS_ISEC } from '../modules/local/bcftools/isec'
include { SNPEFF_SNPEFF } from '../modules/local/snpeff/snpeff'
include { MULTIQC } from '../modules/local/multiqc'

//
// PLUGINS & UTILS: nf-core standard utilities
//
include { paramsSummaryMap } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from 'plugin/nf-core-utils'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { paramsSummaryMultiqc } from '../subworkflows/nf-core/utils_nfcore_pipeline'

//
// SUBWORKFLOWS: Local processing logic
//
include { FASTQ_FASTQC_EXTRACT_CUTADAPT_BBSPLIT_SORTMERNA } from '../subworkflows/local/fastq_fastqc_extract_cutadapt_bbsplit_sortmerna'
include { FASTQ_SAMPLE_SEQKIT_SALMON } from '../subworkflows/local/fastq_sample_seqkit_salmon'
include { BAM_SORT_INDEX_DEDUP_SAMTOOLS_UMITOOLS_FEATURECOUNTS } from '../subworkflows/local/bam_sort_index_dedup_samtools_umitools_featurecounts'
include { BAM_RSEQC } from '../subworkflows/local/bam_rseqc'
include { BAM_DEXSEQ_DEU_ENRICHMENT } from '../subworkflows/local/bam_dexseq_deu_enrichment'
include { BAM_BAMLIST_AVGLEN_RMATS } from '../subworkflows/local/bam_bamlist_avglen_rmats'
include { FASTQ_QUANT_PSEUDOALIGNMENT } from '../subworkflows/local/fastq_quant_pseudoalignment'
include { COUNTS_MATRIXGENERATION } from '../subworkflows/local/counts_matrixgeneration'
include { COUNTS_COEXPRESSION_NETWORK_PPI_ENRICHMENT } from '../subworkflows/local/counts_coexpression_network_ppi_enrichment'
include { COUNTS_DIFFEXPR_QC_ENRICHMENT_VISUALIZATION } from '../subworkflows/local/counts_diffexpr_qc_enrichment_visualization'
include { COUNTS_ISOFORM_SWITCH_ENRICHMENT } from '../subworkflows/local/counts_isoform_switch_enrichment'
include { BAM_DEDUP_RECALIBRATE_PICARD_GATK } from '../subworkflows/local/bam_dedup_recalibrate_picard_gatk'
include { VCF_CALL_FILTER_GATK } from '../subworkflows/local/vcf_call_filter_gatk'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASEQ {
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
    ch_hisat2_index // channel: [ hisat2_index ]
    ch_kallisto_index // channel: [ kallisto_index ]
    ch_salmon_index // channel: [ salmon_index ]
    ch_sortmerna_index // channel: [ sortmerna_index ]
    ch_bbsplit_index // channel: [ bbsplit_index ]
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

    // -------------------------------------------------------------------------
    //  STAGE 0: Initialize and Pre-process
    // -------------------------------------------------------------------------
    ch_input = ch_samplesheet

    // Select appropriate adapter file based on endedness
    ch_adapters = ch_input
        .map { meta_map, fastqs ->
            meta_map.single_end
                ? [[id: "adapters-SE"], file("${projectDir}/assets/adapters-SE.fa")]
                : [[id: "adapters-PE"], file("${projectDir}/assets/adapters-PE.fa")]
        }
        .first()

    // QC, UMI extraction, Trimming, and Contaminant Removal (BBSplit/SortMeRNA)
    FASTQ_FASTQC_EXTRACT_CUTADAPT_BBSPLIT_SORTMERNA(
        ch_input,
        ch_adapters,
        ch_rrna_db_fasta,
        ch_bbsplit_index,
        ch_sortmerna_index,
        params.rrna_db_type,
        params.with_umi,
        params.skip_umi_extract,
        params.skip_fastqc,
        params.skip_trimming,
        params.skip_bbsplit,
        params.skip_sortmerna,
    )
    ch_reads = FASTQ_FASTQC_EXTRACT_CUTADAPT_BBSPLIT_SORTMERNA.out.reads

    // -------------------------------------------------------------------------
    //  STAGE 1: Strandedness Detection (Auto-detection via Salmon)
    // -------------------------------------------------------------------------
    ch_subsampled_reads = channel.empty()
    FASTQ_SAMPLE_SEQKIT_SALMON(
        ch_reads.filter { meta, reads -> meta.lib_type == "auto" },
        ch_gtf_compressed,
        ch_salmon_index,
    )

    ch_subsampled_reads = FASTQ_SAMPLE_SEQKIT_SALMON.out.reads.first()

    // Update metadata with inferred or user-defined strandedness
    ch_reads = FASTQ_SAMPLE_SEQKIT_SALMON.out.strandness
        .join(ch_reads, remainder: true)
        .map { meta, strandness, reads ->
            def updated_meta = meta + [
                lib_type: meta.lib_type in ["forward", "reverse"] ? meta.lib_type : strandness
            ]
            return [updated_meta, reads]
        }

    // Initialize channels for specific alignment outputs
    ch_hisat2_summary = channel.empty()
    ch_bam = channel.empty()
    ch_bai = channel.empty()
    ch_umi_log = channel.empty()
    ch_counts = channel.empty()
    ch_counts_summary = channel.empty()
    ch_genebodycoverage = channel.empty()
    ch_bamstat = channel.empty()
    ch_innerdistance = channel.empty()
    ch_inferexperiment = channel.empty()
    ch_junctionannotation = channel.empty()
    ch_readdistribution = channel.empty()
    ch_readduplication = channel.empty()
    ch_tin = channel.empty()

    // -------------------------------------------------------------------------
    //  STAGE 2: Alignment and Quantification
    // -------------------------------------------------------------------------

    // Case A: Standard Alignment (HISAT2)
    if (params.aligner == "hisat2") {

        HISAT2_ALIGN(
            ch_reads,
            ch_hisat2_index,
            ch_fasta_uncompressed,
        )
        ch_hisat2_summary = HISAT2_ALIGN.out.summary

        // Post-alignment processing: Sort, Index, Dedup, and FeatureCounts
        BAM_SORT_INDEX_DEDUP_SAMTOOLS_UMITOOLS_FEATURECOUNTS(
            HISAT2_ALIGN.out.bam,
            ch_gtf_uncompressed,
            params.with_umi,
            params.analysis_method.split(","),
        )
        ch_bam = BAM_SORT_INDEX_DEDUP_SAMTOOLS_UMITOOLS_FEATURECOUNTS.out.bam
        ch_bai = BAM_SORT_INDEX_DEDUP_SAMTOOLS_UMITOOLS_FEATURECOUNTS.out.bai
        ch_counts = BAM_SORT_INDEX_DEDUP_SAMTOOLS_UMITOOLS_FEATURECOUNTS.out.counts
        ch_counts_summary = BAM_SORT_INDEX_DEDUP_SAMTOOLS_UMITOOLS_FEATURECOUNTS.out.summary
        ch_umi_log = BAM_SORT_INDEX_DEDUP_SAMTOOLS_UMITOOLS_FEATURECOUNTS.out.umi_log

        // RSeQC quality control metrics
        BAM_RSEQC(
            params.rseqc_modules ? params.rseqc_modules.split(",") : channel.empty(),
            ch_bam,
            ch_bed,
            BAM_SORT_INDEX_DEDUP_SAMTOOLS_UMITOOLS_FEATURECOUNTS.out.bai,
        )
        ch_genebodycoverage = BAM_RSEQC.out.genebodycoverage_txt
        ch_readduplication = BAM_RSEQC.out.readduplication_pos
        ch_innerdistance = BAM_RSEQC.out.innerdistance_freq
        ch_bamstat = BAM_RSEQC.out.bamstat
        ch_junctionannotation = BAM_RSEQC.out.junctionannotation_log
        ch_inferexperiment = BAM_RSEQC.out.inferexperiment
        ch_readdistribution = BAM_RSEQC.out.readdistribution
        ch_tin = BAM_RSEQC.out.tin_txt
    }
    else if (params.pseudo_aligner in ["salmon", "kallisto"]) {

        FASTQ_QUANT_PSEUDOALIGNMENT(
            ch_reads,
            ch_kallisto_index,
            ch_salmon_index,
            ch_gtf_compressed,
            params.fragment_length,
            params.fragment_length_sd,
            params.pseudo_aligner,
        )
        ch_counts = FASTQ_QUANT_PSEUDOALIGNMENT.out.quant_dir
        ch_counts_summary = FASTQ_QUANT_PSEUDOALIGNMENT.out.quant_log
    }
    else {
        error(
            """
            ❌ INVALID ALIGNER SELECTION
            Valid combinations: 1. --aligner 'hisat2' | 2. --pseudo_aligner 'salmon' | 3. --pseudo_aligner 'kallisto'
            """
        )
    }

    // -------------------------------------------------------------------------
    //  STAGE 3: Secondary Downstream Analysis (DEG, WGCNA, DIU, AS, DEU)
    // -------------------------------------------------------------------------

    // Count Matrix Generation for Gene-level analysis
    if ("DEG" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",")) {
        COUNTS_MATRIXGENERATION(
            ch_counts.map { meta, counts -> counts }.collect(),
            params.aligner ? "featurecounts" : params.pseudo_aligner,
            params.species,
        )
    }

    // Differential Gene Expression (DEG)
    if ("DEG" in params.analysis_method.split(",")) {
        COUNTS_DIFFEXPR_QC_ENRICHMENT_VISUALIZATION(
            COUNTS_MATRIXGENERATION.out.gene_table_rds,
            params.input,
            params.pvalue_threshold,
            params.logfc_threshold,
            params.ntop_genes,
            params.diffexpr_method,
            params.species,
            params.nes_threshold,
            params.msigdb_categories,
            params.padj_gsea,
            params.ntop_processes,
            params.rank_method,
            params.rsq_threshold,
            params.cluster_method,
            params.enrichment_method ? params.enrichment_method.split(",") : null,
        )
    }

    // Gene Co-expression Network (WGCNA)
    if ("WGCNA" in params.analysis_method.split(",")) {
        COUNTS_COEXPRESSION_NETWORK_PPI_ENRICHMENT(
            COUNTS_MATRIXGENERATION.out.gene_table_rds,
            params.input,
            params.min_gs,
            params.min_mm,
            params.ntop_hubgenes,
            params.ntop_processes,
            params.cor_threshold,
            params.pvalue_threshold,
            params.tomtype,
            params.networktype,
            params.deepsplit,
            params.reassignthreshold,
            params.mergecutheight,
            params.minmodulesize,
            params.sft_r2_threshold,
            params.species,
            params.score_threshold,
            params.enrichment_method ? params.enrichment_method.split(",") : null,
        )
    }

    // Differential Isoform Usage (DIU) & Alternative Splicing (AS - Pseudo)
    if (("DIU" in params.analysis_method.split(",") || "AS" in params.analysis_method.split(",")) && params.pseudo_aligner in ["salmon", "kallisto"]) {
        COUNTS_ISOFORM_SWITCH_ENRICHMENT(
            ch_counts.map { meta, counts -> counts }.collect(),
            params.input,
            ch_gtf_isoform,
            ch_transcriptome,
            params.pseudo_aligner,
            params.analysis_method,
            params.pvalue_threshold,
            params.dif_cutoff,
            params.ntop_isoforms,
            params.species,
            params.ntop_processes,
            params.enrichment_method ? params.enrichment_method.split(",") : null,
        )
    }

    // Alternative Splicing via rMATS (AS - HISAT2 only)
    if ("AS" in params.analysis_method.split(",") && params.aligner == "hisat2") {
        if (ch_subsampled_reads.count() == 0) {
            SEQKIT_SAMPLE(ch_reads.first())
            ch_subsampled_reads = SEQKIT_SAMPLE.out.reads
        }
        BAM_BAMLIST_AVGLEN_RMATS(
            ch_bam,
            ch_subsampled_reads,
            ch_gtf_uncompressed,
            params.statoff,
            params.novelss,
            params.allowclipping,
            params.individualcounts,
            params.event_types,
            params.pvalue_threshold,
            params.delta_psi,
        )
    }

    // Differential Exon Usage (DEU)
    if ("DEU" in params.analysis_method.split(",") && params.aligner == "hisat2") {
        BAM_DEXSEQ_DEU_ENRICHMENT(
            ch_bam,
            ch_gff,
            params.input,
            params.alignment_quality,
            params.pvalue_threshold,
            params.logfc_threshold,
            params.ntop_genes,
            params.min_exonlength,
            params.species,
            params.ntop_processes,
            params.enrichment_method ? params.enrichment_method.split(",") : null,
        )
    }

    // Germline Variant Calling via GATK4 (GVC - HISAT2 only)

    // Initialize channels for specific variant calling outputs
    ch_final_vcf = channel.empty()
    ch_picard_metrics = channel.empty()
    ch_recalibrate_table = channel.empty()

    if ("GVC" in params.analysis_method.split(",") && params.aligner == "hisat2") {
        BAM_DEDUP_RECALIBRATE_PICARD_GATK(
            ch_bam,
            ch_bai,
            ch_fasta_uncompressed,
            ch_fai,
            ch_dict,
            ch_intervals_split,
            ch_knownsites,
            ch_knownsites_tbi,
            params.skip_picard_markduplicates,
            params.skip_baserecalibration,
            params.skip_interval_splitting,
        )
        ch_bam = BAM_DEDUP_RECALIBRATE_PICARD_GATK.out.bam
        ch_bai = BAM_DEDUP_RECALIBRATE_PICARD_GATK.out.bai
        ch_picard_metrics = BAM_DEDUP_RECALIBRATE_PICARD_GATK.out.picard_metrics
        ch_recalibrate_table = BAM_DEDUP_RECALIBRATE_PICARD_GATK.out.recalibrate_table

        VCF_CALL_FILTER_GATK(
            ch_bam,
            ch_bai,
            ch_fasta_uncompressed,
            ch_fai,
            ch_dict,
            ch_intervallist,
            ch_bed,
            ch_dbsnp,
            ch_dbsnp_tbi,
            params.skip_variantfiltration,
        )
        ch_vcf = VCF_CALL_FILTER_GATK.out.vcf
        ch_vcf_tbi = VCF_CALL_FILTER_GATK.out.tbi

        ch_isec_input = ch_vcf
            .join(ch_vcf_tbi)
            .map { meta, vcf, tbi -> [1, vcf, tbi, meta.id] }
            .groupTuple()
            .map { dummy, vcfs, tbis, ids ->
                [[id: ids.join('_')], vcfs, tbis]
            }

        BCFTOOLS_ISEC(
            ch_isec_input,
            [[:], []],
            ch_bed,
        )

        ch_unique_vcf = BCFTOOLS_ISEC.out.results.map { meta, folder ->
            // 1. Identify which condition is NOT control
            def conds = meta.id.split('_')
            def treatment = conds.find { it != 'control' }
            def control = conds.find { it == 'control' }

            // 2. Find the VCF unique to the treatment
            def unique_file = folder
                .listFiles()
                .find {
                    it.name.contains("${treatment}") && it.name.contains('unique') && !it.name.endsWith('.tbi')
                }

            // 3. Create a highly descriptive ID
            def new_id = "${treatment}_absent_in_${control}"

            return [[id: new_id], unique_file]
        }
        ch_vcf_all = ch_vcf.mix(ch_unique_vcf)

        SNPEFF_SNPEFF(
            ch_vcf_all,
            ch_snpeff_genome,
            ch_snpeff_db,
        )
        ch_final_vcf = SNPEFF_SNPEFF.out.vcf
    }

    // -------------------------------------------------------------------------
    //  STAGE 4: MultiQC and Versioning
    // -------------------------------------------------------------------------

    // Software Version Handling
    softwareVersionsToYAML(
        softwareVersions: channel.topic("versions"),
        nextflowVersion: workflow.nextflow.version,
    ).collectFile(
        storeDir: "${params.outdir}/pipeline_info",
        name: 'rnaseq_software_' + 'mqc_' + 'versions.yml',
        sort: true,
        newLine: true,
    ).set { ch_collated_versions }

    ch_multiqc_report = channel.empty()

    if (!params.skip_multiqc) {
        ch_multiqc_config = [[:], file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)]

        // Compile Summary and Description
        summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
        ch_multiqc_files = channel.empty().mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

        ch_multiqc_custom_methods_description = params.multiqc_methods_description
            ? file(params.multiqc_methods_description, checkIfExists: true)
            : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)

        ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
        ch_multiqc_files = ch_multiqc_files.mix(
            ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true)
        )

        // Collate Pre-processing QC
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_EXTRACT_CUTADAPT_BBSPLIT_SORTMERNA.out.fastqc_raw_zip.collect { it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_EXTRACT_CUTADAPT_BBSPLIT_SORTMERNA.out.fastqc_trim_zip.collect { it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_EXTRACT_CUTADAPT_BBSPLIT_SORTMERNA.out.umi_log.collect { it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_EXTRACT_CUTADAPT_BBSPLIT_SORTMERNA.out.trim_json.collect { it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_EXTRACT_CUTADAPT_BBSPLIT_SORTMERNA.out.bbsplit_stats.collect { it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_EXTRACT_CUTADAPT_BBSPLIT_SORTMERNA.out.rrna_log.collect { it[1] }.ifEmpty([]))

        // Collate Alignment/Quant QC
        ch_multiqc_files = ch_multiqc_files.mix(ch_hisat2_summary.collect { it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_umi_log.collect { it[1] }.ifEmpty([]))
        if (params.pseudo_aligner == "salmon") {
            ch_multiqc_files = ch_multiqc_files.mix(ch_counts.collect { it[1] }.ifEmpty([]))
        }
        ch_multiqc_files = ch_multiqc_files.mix(ch_counts_summary.collect { it[1] }.ifEmpty([]))

        // Collate RSeQC metrics
        ch_multiqc_files = ch_multiqc_files.mix(ch_genebodycoverage.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_bamstat.collect { it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_innerdistance.collect { it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_junctionannotation.collect { it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_readduplication.collect { it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_inferexperiment.collect { it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_readdistribution.collect { it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_tin.collect { it[1] }.ifEmpty([]))

        // Collate vairant calling metrics
        ch_multiqc_files = ch_multiqc_files.mix(ch_picard_metrics.collect { it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_recalibrate_table.collect { it[1] }.ifEmpty([]))

        MULTIQC(
            ch_multiqc_files.collect().map { [[:], it] },
            ch_multiqc_config,
        )
        ch_multiqc_report = MULTIQC.out.report
    }

    emit:
    multiqc_report = ch_multiqc_report // channel: /path/to/multiqc_report.html
}
