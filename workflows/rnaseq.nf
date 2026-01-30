/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { HISAT2_ALIGN } from '../modules/local/hisat2/align'
include { GXF2BED } from '../modules/local/gxf2bed'
include { SEQKIT_SAMPLE } from '../modules/local/seqkit/sample'
include { MULTIQC } from '../modules/local/multiqc'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { paramsSummaryMap } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_rnaseq_pipeline'

//
// SUBWORKFLOW: Loaded from subworkflows/local/
//
include { FASTQ_FASTQC_EXTRACT_CUTADAPT_BBSPLIT_SORTMERNA } from '../subworkflows/local/fastq_fastqc_extract_cutadapt_bbsplit_sortmerna'
include { FASTQ_SAMPLE_SEQKIT_SALMON } from '../subworkflows/local/fastq_sample_seqkit_salmon'
include { GTF_EXTRACTSPLICESITES_EXTRACTEXONS_BUILD_HISAT2 } from '../subworkflows/local/gtf_extractsplicesites_extractexons_build_hisat2'
include { BAM_SORT_INDEX_DEDUP_SAMTOOLS_UMITOOLS_FEATURECOUNTS } from '../subworkflows/local/bam_sort_index_dedup_samtools_umitools_featurecounts'
include { BAM_RSEQC } from '../subworkflows/local/bam_rseqc'
include { BAM_DEXSEQ_DEU_ENRICHMENT } from '../subworkflows/local/bam_dexseq_deu_enrichment'
include { BAM_BAMLIST_AVGLEN_RMATS } from '../subworkflows/local/bam_bamlist_avglen_rmats'
include { FASTQ_INDEX_QUANT_PSEUDOALIGNMENT } from '../subworkflows/local/fastq_index_quant_pseudoalignment'
include { COUNTS_MATRIXGENERATION } from '../subworkflows/local/counts_matrixgeneration'
include { COUNTS_COEXPRESSION_NETWORK_PPI_ENRICHMENT } from '../subworkflows/local/counts_coexpression_network_ppi_enrichment'
include { COUNTS_DIFFEXPR_QC_ENRICHMENT_VISUALIZATION } from '../subworkflows/local/counts_diffexpr_qc_enrichment_visualization'
include {
    COUNTS_ISOFORM_SWITCH_ENRICHMENT as COUNTS_DIU_ENRICHMENT ;
    COUNTS_ISOFORM_SWITCH_ENRICHMENT as COUNTS_AS_ENRICHMENT
} from '../subworkflows/local/counts_isoform_switch_enrichment'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASEQ {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_rrna_db
    ch_genome_fasta
    ch_contaminant_fasta
    ch_transcript_fasta
    ch_gtf

    main:

    ch_versions = channel.empty()

    ch_input = channel.fromPath(ch_samplesheet)
        .splitCsv(header: true)
        .map { row ->
            // row is a Map, not individual params
            def meta = row.sample_id
            def fastq_1 = row.fastq_1
            def fastq_2 = row.fastq_2
            def condition = row.condition
            def lib_type = row.lib_type
            def sequencer = row.sequencer

            // Create the complete meta map
            def meta_map = [
                id: meta,
                condition: condition,
                lib_type: lib_type,
                sequencer: sequencer,
            ]

            // Create fastqs list - 1 or 2 files
            def fastqs
            if (!fastq_2) {
                // Single-end
                meta_map.single_end = true
                fastqs = [fastq_1]
            }
            else {
                // Paired-end
                meta_map.single_end = false
                fastqs = [fastq_1, fastq_2]
            }

            // Return exactly 2 elements: [meta_map, fastqs]
            return [meta_map, fastqs]
        }

    ch_adapters = ch_input
        .map { meta_map, fastqs ->
            def adapter_file = meta_map.single_end
                ? file("${projectDir}/assets/adapters-SE.fa")
                : file("${projectDir}/assets/adapters-PE.fa")
            return adapter_file
        }
        .first()

    FASTQ_FASTQC_EXTRACT_CUTADAPT_BBSPLIT_SORTMERNA(
        ch_input,
        ch_adapters,
        ch_rrna_db,
        ch_genome_fasta,
        ch_contaminant_fasta,
        params.rrna_db_type,
        params.lib_kit,
        params.with_umi,
        params.skip_umi_extract,
        params.skip_fastqc,
        params.skip_trimming,
        params.skip_bbsplit,
        params.skip_sortmerna,
    )
    ch_reads = FASTQ_FASTQC_EXTRACT_CUTADAPT_BBSPLIT_SORTMERNA.out.reads
    ch_versions = ch_versions.mix(FASTQ_FASTQC_EXTRACT_CUTADAPT_BBSPLIT_SORTMERNA.out.versions)

    FASTQ_SAMPLE_SEQKIT_SALMON(
        ch_reads.filter { meta, reads ->
            meta.lib_type == "auto"
        },
        ch_transcript_fasta,
        ch_gtf,
    )
    ch_versions = ch_versions.mix(FASTQ_SAMPLE_SEQKIT_SALMON.out.versions)

    ch_reads = FASTQ_SAMPLE_SEQKIT_SALMON.out.strandness
        .join(
            ch_reads,
            remainder: true
        )
        .map { meta, strandness, reads ->
            def updated_meta = meta + [
                lib_type: meta.lib_type in ["forward", "reverse"] ? meta.lib_type : strandness
            ]
            [updated_meta, reads]
        }

    ch_hisat2_summary = channel.empty()
    ch_bam = channel.empty()
    ch_umi_log = channel.empty()
    ch_counts = channel.empty()
    ch_counts_summary = channel.empty()

    ch_genebodycoverage = channel.empty()
    ch_bamstat = channel.empty()
    ch_inferexperiment = channel.empty()
    ch_readdistribution = channel.empty()
    ch_tin = channel.empty()

    if (params.aligner == "hisat2" || "DEU" in params.analysis_method.split(",")) {
        GTF_EXTRACTSPLICESITES_EXTRACTEXONS_BUILD_HISAT2(
            ch_genome_fasta,
            ch_gtf,
        )
        ch_versions = ch_versions.mix(GTF_EXTRACTSPLICESITES_EXTRACTEXONS_BUILD_HISAT2.out.versions)

        HISAT2_ALIGN(
            ch_reads,
            GTF_EXTRACTSPLICESITES_EXTRACTEXONS_BUILD_HISAT2.out.index,
            ch_genome_fasta,
        )
        ch_hisat2_summary = HISAT2_ALIGN.out.summary
        ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions.first())

        BAM_SORT_INDEX_DEDUP_SAMTOOLS_UMITOOLS_FEATURECOUNTS(
            HISAT2_ALIGN.out.bam,
            ch_gtf,
            params.with_umi,
        )
        ch_bam = BAM_SORT_INDEX_DEDUP_SAMTOOLS_UMITOOLS_FEATURECOUNTS.out.bam
        ch_counts = BAM_SORT_INDEX_DEDUP_SAMTOOLS_UMITOOLS_FEATURECOUNTS.out.counts
        ch_counts_summary = BAM_SORT_INDEX_DEDUP_SAMTOOLS_UMITOOLS_FEATURECOUNTS.out.summary
        ch_umi_log = BAM_SORT_INDEX_DEDUP_SAMTOOLS_UMITOOLS_FEATURECOUNTS.out.umi_log
        ch_versions = ch_versions.mix(BAM_SORT_INDEX_DEDUP_SAMTOOLS_UMITOOLS_FEATURECOUNTS.out.versions)

        GXF2BED(
            ch_gtf
        )
        ch_versions = ch_versions.mix(GXF2BED.out.versions)

        BAM_RSEQC(
            params.rseqc_modules.split(","),
            ch_bam,
            GXF2BED.out.bed,
            BAM_SORT_INDEX_DEDUP_SAMTOOLS_UMITOOLS_FEATURECOUNTS.out.bai,
        )
        ch_genebodycoverage = BAM_RSEQC.out.genebodycoverage_txt
        ch_bamstat = BAM_RSEQC.out.bamstat
        ch_inferexperiment = BAM_RSEQC.out.inferexperiment
        ch_readdistribution = BAM_RSEQC.out.readdistribution
        ch_tin = BAM_RSEQC.out.tin_txt
        ch_versions = ch_versions.mix(BAM_RSEQC.out.versions)
    }
    else if (params.pseudo_aligner in ["salmon", "kallisto"] || "DIU" in params.analysis_method.split(",")) {

        FASTQ_INDEX_QUANT_PSEUDOALIGNMENT(
            ch_reads,
            ch_transcript_fasta,
            ch_gtf,
            params.bootstrap_count,
            params.fragment_length,
            params.fragment_length_sd,
            params.pseudo_aligner,
        )
        ch_counts = FASTQ_INDEX_QUANT_PSEUDOALIGNMENT.out.quant_dir
        ch_counts_summary = FASTQ_INDEX_QUANT_PSEUDOALIGNMENT.out.quant_log
        ch_versions = ch_versions.mix(FASTQ_INDEX_QUANT_PSEUDOALIGNMENT.out.versions)
    }
    else {
        error(
            """
            ❌ INVALID ALIGNER SELECTION

            You specified: --aligner '${params.aligner}' --pseudo_aligner '${params.pseudo_aligner}'

            Valid combinations are:
            1. --aligner 'hisat2'         (HISAT2)
            2. --pseudo_aligner 'salmon'         (Salmon)
            3. --pseudo_aligner 'kallisto'       (Kallisto)

            Please choose one of the above valid options.
            """
        )
    }

    COUNTS_MATRIXGENERATION(
        ch_counts,
        params.quant_method,
        params.species,
    )
    ch_versions = ch_versions.mix(COUNTS_MATRIXGENERATION.out.versions)
    if ("DEG" in params.analysis_method.split(",")) {
        COUNTS_DIFFEXPR_QC_ENRICHMENT_VISUALIZATION(
            COUNTS_MATRIXGENERATION.out.gene_table_rds,
            ch_samplesheet,
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
            params.enrichment_method.split(","),
        )
        ch_versions = ch_versions.mix(COUNTS_DIFFEXPR_QC_ENRICHMENT_VISUALIZATION.out.versions)
    }
    else if ("WGCNA" in params.analysis_method.split(",")) {
        COUNTS_COEXPRESSION_NETWORK_PPI_ENRICHMENT(
            COUNTS_MATRIXGENERATION.out.gene_table_rds,
            ch_samplesheet,
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
            params.species,
            params.score_threshold,
            params.enrichment_method.split(","),
        )
        ch_versions = ch_versions.mix(COUNTS_COEXPRESSION_NETWORK_PPI_ENRICHMENT.out.versions)
    }
    else if ("DIU" in params.analysis_method.split(",")) {
        COUNTS_DIU_ENRICHMENT(
            ch_counts,
            ch_samplesheet,
            ch_gtf,
            ch_transcript_fasta,
            params.pseudo_aligner,
            "DIU",
            params.pvalue_threshold,
            params.dif_cutoff,
            params.ntop_isoforms,
            params.species,
            params.ntop_processes,
            params.enrichment_method.split(","),
        )
        ch_versions = ch_versions.mix(COUNTS_DIU_ENRICHMENT.out.versions)
    }
    else if ("AS" in params.analysis_method.split(",")) {

        SEQKIT_SAMPLE(
            ch_reads.first()
        )
        ch_versions = ch_versions.mix(SEQKIT_SAMPLE.out.versions)

        BAM_BAMLIST_AVGLEN_RMATS(
            ch_bam,
            SEQKIT_SAMPLE.out.reads,
            ch_gtf,
            params.statoff,
            params.novelss,
            params.allowclipping,
            params.individualcounts,
            params.event_types,
            params.pvalue_threshold,
            params.delta_psi,
        )
        ch_versions = ch_versions.mix(BAM_BAMLIST_AVGLEN_RMATS.out.versions)

        COUNTS_AS_ENRICHMENT(
            ch_counts,
            ch_samplesheet,
            ch_gtf,
            ch_transcript_fasta,
            params.pseudo_aligner,
            "AS",
            params.pvalue_threshold,
            params.dif_cutoff,
            params.ntop_isoforms,
            params.species,
            params.ntop_processes,
            params.enrichment_method.split(","),
        )
        ch_versions = ch_versions.mix(COUNTS_AS_ENRICHMENT.out.versions)
    }
    else {
        error(
            """
            ❌ INVALID ANALYSIS METHOD(S)

            Valid analysis methods are:
              • DEG   - Differential Expression (limma/DESeq2/edgeR)
              • WGCNA - Co-expression Network Analysis
              • DIU   - Differential Isoform Usage
              • AS    - Alternative Splicing (rMATS)

            You requested: ${params.analysis_method}

            Example usage:
              --analysis_method DEG,WGCNA   # Run DE and co-expression
              --analysis_method DIU,AS      # Run isoform and splicing
              --analysis_method DEG         # Run only differential expression
            """
        )
    }

    //
    // Collate and save software versions
    //
    def topic_versions = Channel
        .topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}"]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'rnaseq_software_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_report = channel.empty()

    if (!params.skip_multiqc) {
        ch_multiqc_files = channel.empty()

        // Mix tool outputs - extract just the file paths
        ch_multiqc_files = ch_multiqc_files.mix(ch_hisat2_summary.map { meta, summary -> summary }.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_umi_log.map { meta, log -> log }.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_counts_summary.map { meta, summary -> summary }.collect().ifEmpty([]))

        // Mix RSeQC outputs
        ch_multiqc_files = ch_multiqc_files.mix(ch_genebodycoverage.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_bamstat.map { meta, bamstat -> bamstat }.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_inferexperiment.map { meta, infer -> infer }.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_readdistribution.map { meta, dist -> dist }.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_tin.map { meta, tin -> tin }.collect().ifEmpty([]))

        MULTIQC(
            ch_multiqc_files.collect()
        )
        ch_multiqc_report = MULTIQC.out.report
    }

    emit:
    multiqc_report = ch_multiqc_report // channel: /path/to/multiqc_report.html
    versions = ch_versions // channel: [ path(versions.yml) ]
}
