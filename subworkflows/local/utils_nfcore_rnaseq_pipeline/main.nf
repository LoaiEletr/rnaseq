//
// Subworkflow with functionality specific to the LoaiEletr/rnaseq pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap } from 'plugin/nf-schema'
include { samplesheetToList } from 'plugin/nf-schema'
include { paramsHelp } from 'plugin/nf-schema'
include { completionEmail } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {
    take:
    version // boolean: Display version and exit
    validate_params // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir //  string: The output directory where the results will be saved
    input //  string: Path to input samplesheet
    help // boolean: Display help message and exit
    help_full // boolean: Show the full help message
    show_hidden // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1,
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  Loai3tr/comprehensive-rnaseq-pipeline ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/', '')}" }.join("\n")}${workflow.manifest.doi ? "\n" : ""}

* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/Loai3tr/comprehensive-rnaseq-pipeline/blob/master/CITATIONS.md
"""
    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --species human --genome GRCh38 --contaminant_species mouse --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN(
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE(
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from YOUR CUSTOM input samplesheet format
    // MATCHING YOUR WORKFLOW'S EXPECTATIONS EXACTLY
    //
    ch_samplesheet = channel.fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map { meta, fastq_1, fastq_2 ->
            if (!fastq_2) {
                return [meta + [single_end: true], [fastq_1]]
            }
            else {
                return [meta + [single_end: false], [fastq_1, fastq_2]]
            }
        }
    ch_samplesheet
        .map { meta, fastq_list -> meta }
        .collect()
        .subscribe { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }

    emit:
    samplesheet = ch_samplesheet // channel: sample fastqs parsed from --input
    versions = ch_versions // channel: [ versions.yml ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {
    take:
    email //  string: email address
    email_on_fail //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url //  string: hook URL for notifications
    multiqc_report //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error("Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters for RNA-seq pipeline
//
def validateInputParameters() {
    genomeExistsError()
    contaminantExistsError()

    // Validate aligner selection: exactly one must be provided and valid
    def has_aligner = ["hisat2"].contains(params.aligner)
    def has_pseudo = ["salmon", "kallisto"].contains(params.pseudo_aligner)

    if (has_aligner && has_pseudo) {
        error("❌ CONFLICTING ALIGNERS: Please specify either --aligner OR --pseudo_aligner, not both.")
    }

    if (!has_aligner && !has_pseudo) {
        error(
            """
    ❌ INVALID OR MISSING ALIGNER SELECTION

    You must specify exactly one valid method:
      --aligner 'hisat2'
      OR
      --pseudo_aligner 'salmon'
      --pseudo_aligner 'kallisto'

    Current selection: --aligner '${params.aligner}' --pseudo_aligner '${params.pseudo_aligner}'
    """
        )
    }

    // Validate analysis methods
    def valid_methods = ["DEG", "WGCNA", "DIU", "AS", "DEU"]
    def requested_methods = params.analysis_method ? params.analysis_method.split(",").collect { it.trim() } : []
    def invalid_methods = requested_methods.findAll { !valid_methods.contains(it) }

    if (invalid_methods) {
        error(
            """
        ❌ INVALID ANALYSIS METHOD(S): ${invalid_methods.join(", ")}

        Valid analysis methods are:
          • DEG   - Differential Expression (limma/DESeq2/maSigPro)
          • WGCNA - Co-expression Network Analysis
          • DIU   - Differential Isoform Usage (ISoformSwitchAnalyzeR)
          • DEU   - Differential Exon Usage (DEXSeq)
          • AS    - Alternative Splicing (rMATS, ISoformSwitchAnalyzeR)

        Example usage:
          --analysis_method DEG,WGCNA   # Run DE and co-expression
          --analysis_method DIU,AS      # Run isoform and splicing
        """
        )
    }

    // Validate RSeQC modules
    def valid_rseqc = ['bam_stat', 'genebody_coverage', 'infer_experiment', 'inner_distance', 'junction_annotation', 'read_distribution', 'read_duplication', 'tin']
    def requestedRseqc = params.rseqc_modules ? params.rseqc_modules.split(',').collect { it.trim() } : []
    def invalid_rseqc = requestedRseqc.findAll { !valid_rseqc.contains(it) }

    if (invalid_rseqc) {
        error("❌ INVALID RSeQC MODULE(S): ${invalid_rseqc.join(', ')}\nValid modules are: ${valid_rseqc.join(', ')}")
    }

    // Validate MSigDB categories
    def valid_msigdb = ['H', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
    def requested_msigdb = params.msigdb_categories ? params.msigdb_categories.split(',').collect { it.trim() } : []
    def invalid_msigdb = requested_msigdb.findAll { !valid_msigdb.contains(it) }

    if (invalid_msigdb) {
        error(
            """
        ❌ INVALID MSigDB CATEGORY: ${invalid_msigdb.join(', ')}

        Valid categories are:
          ${valid_msigdb.join(', ')}
        """
        )
    }

    // Validate WGCNA parameters
    def valid_tom_types = ["none", "signed", "signed 2", "unsigned", "unsigned 2", "signed Nowick", "signed Nowick 2"]
    def valid_network_types = ["signed", "unsigned", "signed hybrid"]

    if ("WGCNA" in requested_methods) {
        if (!valid_tom_types.contains(params.tomtype)) {
            error(
                """
            ❌ INVALID WGCNA TOM TYPE: ${params.tomtype}

            Valid types are:
              ${valid_tom_types.join(', ')}
            """
            )
        }
        if (!valid_network_types.contains(params.networktype)) {
            error(
                """
            ❌ INVALID WGCNA NETWORK TYPE: ${params.networktype}

            Valid types are:
              ${valid_network_types.join(', ')}
            """
            )
        }
    }

    // Validate Alternative Splicing (AS) event types
    def valid_events = ["SE", "RI", "MXE", "A5SS", "A3SS"]
    def requested_events = params.event_types ? params.event_types.split(',').collect { it.trim() } : []
    def invalid_events = requested_events.findAll { !valid_events.contains(it) }

    if ("AS" in requested_methods && invalid_events) {
        error(
            """
        ❌ INVALID AS EVENT TYPE(S): ${invalid_events.join(', ')}

        Valid event types are:
          ${valid_events.join(', ')}
        """
        )
    }

    // Validate Differential Expression and Clustering methods
    def valid_diff_methods = ["limma", "deseq2", "masigpro"]
    def valid_cluster_methods = ["hclust", "Mclust", "kmeans"]
    def valid_rank_methods = ["t_stat", "logfc", "signed_significance"]

    if (!valid_diff_methods.contains(params.diffexpr_method)) {
        error("❌ INVALID DIFFEXPR METHOD: ${params.diffexpr_method}\nValid: ${valid_diff_methods.join(', ')}")
    }
    if (!valid_rank_methods.contains(params.rank_method)) {
        error("❌ INVALID RANK METHOD: ${params.rank_method}\nValid: ${valid_rank_methods.join(', ')}")
    }
    if (!valid_cluster_methods.contains(params.cluster_method)) {
        error("❌ INVALID CLUSTER METHOD: ${params.cluster_method}\nValid: ${valid_cluster_methods.join(', ')}")
    }

    // Validate Library and rRNA parameters
    def valid_kits = ["quantseq", "corall", "takara"]
    def valid_rrnas = ["fast", "default", "sensitive", "sensitive_rfam"]

    if (!valid_kits.contains(params.lib_kit)) {
        error("❌ INVALID LIB KIT: ${params.lib_kit}\nOptions: ${valid_kits.join(', ')}")
    }
    if (!valid_rrnas.contains(params.rrna_db_type)) {
        error("❌ INVALID rRNA DB TYPE: ${params.rrna_db_type}\nOptions: ${valid_rrnas.join(', ')}")
    }

    // Validate Species
    def valid_species = ["human", "mouse", "rat", "yeast", "fruitfly", "zebrafish", "worm", "arabidopsis", "chicken", "cow", "pig", "dog", "monkey"]
    if (!valid_species.contains(params.species)) {
        error(
            """
        ❌ INVALID SPECIES: ${params.species}

        Supported species are:
          ${valid_species.join(', ')}
        """
        )
    }
    else if (params.contaminant_species) {
        if (!valid_species.contains(params.contaminant_species)) {
            error(
                """
        ❌ INVALID SPECIES: ${params.contaminant_species}

        Supported species are:
          ${valid_species.join(', ')}
        """
            )
        }
    }

    // Validate Enrichment methods
    def valid_enrichment = ["GO", "KEGG", "GSEA"]
    def requested_enrichment = params.enrichment_method ? params.enrichment_method.split(",").collect { it.trim() } : []
    def invalid_enrichment = requested_enrichment.findAll { !valid_enrichment.contains(it) }

    if (invalid_enrichment) {
        error(
            """
        ❌ INVALID ENRICHMENT METHOD(S): ${invalid_enrichment.join(", ")}

        Valid enrichment methods are:
          • GO   - Gene Ontology Enrichment
          • KEGG - Kyoto Encyclopedia of Genes and Genomes
          • GSEA - Gene Set Enrichment Analysis

        Example usage:
          --enrichment_method GO,KEGG
        """
        )
    }

    // Validate combination of aligner and analysis methods
    if (params.aligner == "hisat2" && "DIU" in requested_methods) {
        error("❌  DIU analysis requires pseudo-aligner (salmon/kallisto), but HISAT2 is selected. Consider using --pseudo_aligner instead.")
    }
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(metas) {
    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect { meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("❌  INVALID INPUT: All samples must be of the same datatype (Single-End or Paired-End). Mixed endedness detected.")
    }

    // Validate that exactly 2 conditions are present
    def conditions = metas.collect { it.condition }.unique()
    if (conditions.size() != 2) {
        error("❌  INSUFFICIENT CONDITIONS: Comparative analysis requires at least 2 distinct conditions. Found: ${conditions}.")
    }

    // Check that each condition has sufficient replicates for robust statistical analysis
    def invalid_groups = metas.countBy { it.condition }.findAll { condition, count -> count < 3 }
    if (invalid_groups) {
        error(
            "❌ INSUFFICIENT REPLICATES: The following groups have less than 3 samples:\n" + invalid_groups.collect { condition, count -> "  • ${condition}: ${count} sample(s)" }.join("\n") + "\nStatistical analysis requires at least 3 replicates per condition for reliable results."
        )
    }

    // Check lib_type consistency
    def lib_types = metas.collect { it.lib_type }.unique()
    if (lib_types.size() > 1) {
        error("❌  INCONSISTENT LIBRARY TYPES: Detected multiple types: ${lib_types}. All samples in a single run must use the same orientation.")
    }

    // Validate lib_type values
    def valid_lib_types = ["forward", "reverse", "auto"]
    def invalid_lib_types = lib_types.findAll { !valid_lib_types.contains(it) }
    if (invalid_lib_types) {
        error("❌ INVALID LIBRARY STRANDEDNESS: '${invalid_lib_types.join(', ')}' is not supported. Please use: forward, reverse, or auto.")
    }

    // Check sequencer information
    def valid_sequencers = ["HiSeq", "MiSeq", "NovaSeq", "NextSeq"]
    def sequencers = metas.collect { it.sequencer }.unique()

    // Check for unsupported platforms
    def invalid_sequencers = sequencers.findAll { !valid_sequencers.contains(it) }
    if (invalid_sequencers) {
        error("❌  UNSUPPORTED SEQUENCER: '${invalid_sequencers.join(', ')}' is not in the supported list. Valid options: ${valid_sequencers.join(', ')}.")
    }

    // Check for mixed platforms
    if (sequencers.size() > 1) {
        error("❌  MIXED SEQUENCING PLATFORMS: Detected multiple sequencers: ${sequencers}. To avoid significant batch effects, please process data from different sequencers in separate runs.")
    }
}
//
// Get attribute from genome config file e.g. fasta, gtf
//
def getGenomeAttribute(attribute) {
    if (params.species && params.genome && params.genomes.containsKey(params.species) && params.genomes[params.species].containsKey(params.genome)) {

        if (params.genomes[params.species][params.genome].containsKey(attribute)) {
            return params.genomes[params.species][params.genome][attribute]
        }
    }
    log.warn("Genome attribute '${attribute}' not found for species=${params.species}, genome=${params.genome}")
    return null
}

//
// Get contaminant genome from contaminant config
//
def getContaminantGenome() {
    if (params.contaminant_species && params.contaminants && params.contaminants.containsKey(params.contaminant_species)) {
        return params.contaminants[params.contaminant_species].genome
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.species && params.genome) {
        // Check if species exists in genomes config
        if (!params.genomes.containsKey(params.species)) {
            def error_string = """
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Species '${params.species}' not found in genome config.
  Available species: ${params.genomes.keySet().join(", ")}
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
            error(error_string)
        }

        // Check if genome exists for that species
        if (params.genomes.containsKey(params.species) && !params.genomes[params.species].containsKey(params.genome)) {

            def available_genomes = params.genomes[params.species].keySet().join(", ")
            def error_string = """
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Genome '${params.genome}' not found for species '${params.species}'.

  Available genomes for ${params.species}:
  ${available_genomes}

  Usage: --species ${params.species} --genome <one_of_the_above>
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
            error(error_string)
        }
    }
}

//
// Exit pipeline if incorrect --contaminant_species provided
//
def contaminantExistsError() {
    if (params.contaminant_species) {
        // Check if contaminant species exists in config
        if (!params.contaminants.containsKey(params.contaminant_species)) {
            def available_contaminants = params.contaminants.keySet().join(", ")
            def error_string = """
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Contaminant species '${params.contaminant_species}' not found in config.

  Available contaminant species:
  ${available_contaminants}

  Usage: --contaminant_species <one_of_the_above>
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
            error(error_string)
        }
    }
}
//
// Generate methods description for MultiQC - RNA-seq specific
//
def toolCitationText() {
    def citation_text = [
        "Tools used in the workflow included:",
        !params.skip_fastqc ? "FastQC (<a href=\"https://www.bioinformatics.babraham.ac.uk/projects/fastqc/\" target=\"_blank\">Andrews, 2010</a>)," : "",
        !params.skip_multiqc ? "MultiQC (<a href=\"https://doi.org/10.1093/bioinformatics/btw354\" target=\"_blank\">Ewels <em>et al.</em>, 2016</a>)," : "",
        "SeqKit (<a href=\"https://doi.org/10.1371/journal.pone.0163962\" target=\"_blank\">Shen <em>et al.</em>, 2016</a>),",
        params.aligner == "hisat2" ? "HISAT2 (<a href=\"https://doi.org/10.1038/s41587-019-0201-4\" target=\"_blank\">Kim <em>et al.</em>, 2019</a>)," : "",
        params.pseudo_aligner == "salmon" ? "Salmon (<a href=\"https://doi.org/10.1038/nmeth.4197\" target=\"_blank\">Patro <em>et al.</em>, 2017</a>)," : "",
        params.pseudo_aligner == "kallisto" ? "Kallisto (<a href=\"https://doi.org/10.1038/nbt.3519\" target=\"_blank\">Bray <em>et al.</em>, 2016</a>)," : "",
        params.with_umi ? "UMI-tools (<a href=\"https://doi.org/10.1101/gr.209601.116\" target=\"_blank\">Smith <em>et al.</em>, 2017</a>)," : "",
        !params.skip_trimming ? "Cutadapt (<a href=\"https://doi.org/10.14806/ej.17.1.200\" target=\"_blank\">Martin, 2011</a>)," : "",
        params.aligner == "hisat2" ? "SAMtools (<a href=\"https://doi.org/10.1093/bioinformatics/btp352\" target=\"_blank\">Li <em>et al.</em>, 2009</a>)," : "",
        params.aligner == "hisat2" && "DEG" in params.analysis_method.split(",") ? "featureCounts (<a href=\"https://doi.org/10.1093/bioinformatics/btt656\" target=\"_blank\">Liao <em>et al.</em>, 2014</a>)," : "",
        "AS" in params.analysis_method.split(",") && params.aligner == "hisat2" ? "rMATS-turbo (<a href=\"https://doi.org/10.1038/s41596-023-00944-2\" target=\"_blank\">Wang <em>et al.</em>, 2024</a>)," : "",
        "AS" in params.analysis_method.split(",") && params.aligner == "hisat2" ? "rMATS2SashimiPlot (<a href=\"https://doi.org/10.5281/ZENODO.10008656\" target=\"_blank\">Shieh <em>et al.</em>, 2023</a>)," : "",
        params.aligner == "hisat2" && (params.rseqc_modules ? params.rseqc_modules.split(",").any { it in ["bam_stat", "genebody_coverage", "infer_experiment", "inner_distance", "junction_annotation", "read_distribution", "read_duplication", "tin"] } : null) ? "RSeQC (<a href=\"https://doi.org/10.1093/bioinformatics/bts356\" target=\"_blank\">Wang <em>et al.</em>, 2012</a>)," : "",
        "DEU" in params.analysis_method.split(",") ? "DEXSeq (<a href=\"https://doi.org/10.1101/gr.133744.111\" target=\"_blank\">Anders <em>et al.</em>, 2012</a>)," : "",
        "DEG" in params.analysis_method.split(",") && (params.diffexpr_method == "limma" || params.diffexpr_method == "masigpro") ? "edgeR (<a href=\"https://doi.org/10.1093/nar/gkaf018\" target=\"_blank\">Chen <em>et al.</em>, 2025</a>)," : "",
        "DEG" in params.analysis_method.split(",") && params.diffexpr_method == "limma" ? "limma (<a href=\"https://doi.org/10.1093/nar/gkv007\" target=\"_blank\">Ritchie <em>et al.</em>, 2015</a>)," : "",
        ("DEG" in params.analysis_method.split(",") && params.diffexpr_method == "deseq2") || "WGCNA" in params.analysis_method.split(",") ? "DESeq2 (<a href=\"https://doi.org/10.1186/s13059-014-0550-8\" target=\"_blank\">Love <em>et al.</em>, 2014</a>)," : "",
        "DEG" in params.analysis_method.split(",") && params.diffexpr_method == "masigpro" ? "reshape (<a href=\"https://doi.org/10.18637/jss.v021.i12\" target=\"_blank\">Wickham, 2007</a>)," : "",
        "DEG" in params.analysis_method.split(",") && params.diffexpr_method == "masigpro" ? "patchwork (<a href=\"https://doi.org/10.32614/cran.package.patchwork\" target=\"_blank\">Pedersen, 2019</a>)," : "",
        "DEG" in params.analysis_method.split(",") && params.diffexpr_method == "masigpro" ? "mclust (<a href=\"https://doi.org/10.1201/9781003277965\" target=\"_blank\">Scrucca <em>et al.</em>, 2023</a>)," : "",
        "DEG" in params.analysis_method.split(",") && params.pseudo_aligner in ["kallisto", "salmon"] ? "tximport (<a href=\"https://doi.org/10.12688/f1000research.7563.2\" target=\"_blank\">Soneson <em>et al.</em>, 2015</a>)," : "",
        "DEG" in params.analysis_method.split(",") && params.pseudo_aligner in ["kallisto", "salmon"] ? "rhdf5 (<a href=\"https://bioconductor.org/packages/rhdf5\" target=\"_blank\">Fischer <em>et al.</em>, 2025</a>)," : "",
        "DEG" in params.analysis_method.split(",") && params.pseudo_aligner in ["kallisto", "salmon"] ? "jsonlite (<a href=\"https://doi.org/10.48550/ARXIV.1403.2805\" target=\"_blank\">Ooms, 2014</a>)," : "",
        "DEG" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || (("DIU" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null)) ? "biomaRt (<a href=\"https://doi.org/10.1038/nprot.2009.97\" target=\"_blank\">Durinck <em>et al.</em>, 2009</a>)," : "",
        "DEG" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") ? "pheatmap (<a href=\"https://github.com/raivokolde/pheatmap\" target=\"_blank\">Kolde, 2025</a>)," : "",
        "DEG" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") ? "EnhancedVolcano (<a href=\"https://github.com/kevinblighe/EnhancedVolcano\" target=\"_blank\">Blighe <em>et al.</em>, 2019</a>)," : "",
        "DEG" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") ? "RColorBrewer (<a href=\"https://doi.org/10.32614/cran.package.rcolorbrewer\" target=\"_blank\">Neuwirth, 2002</a>)," : "",
        "DEG" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || ("AS" in params.analysis_method.split(",") && params.pseudo_aligner in ["salmon", "kallisto"]) ? "ggplot2 (<a href=\"https://ggplot2.tidyverse.org\" target=\"_blank\">Wickham, 2016</a>)," : "",
        "DEG" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || (("DIU" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null)) ? "tibble (<a href=\"https://doi.org/10.32614/cran.package.tibble\" target=\"_blank\">Müller &amp; Wickham, 2016</a>)," : "",
        "DEG" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") ? "tidyr (<a href=\"https://doi.org/10.32614/cran.package.tidyr\" target=\"_blank\">Wickham <em>et al.</em>, 2014</a>)," : "",
        "DEG" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || (("DIU" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null)) ? "dplyr (<a href=\"https://doi.org/10.32614/cran.package.dplyr\" target=\"_blank\">Wickham <em>et al.</em>, 2014</a>)," : "",
        "DEG" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") ? "cowplot (<a href=\"https://doi.org/10.32614/cran.package.cowplot\" target=\"_blank\">Wilke, 2025</a>)," : "",
        "WGCNA" in params.analysis_method.split(",") ? "WGCNA (<a href=\"https://doi.org/10.1186/1471-2105-9-559\" target=\"_blank\">Langfelder &amp; Horvath, 2008</a>)," : "",
        "WGCNA" in params.analysis_method.split(",") ? "STRINGdb (<a href=\"https://doi.org/10.1093/nar/gkac1000\" target=\"_blank\">Szklarczyk <em>et al.</em>, 2023</a>)," : "",
        ("DIU" in params.analysis_method.split(",") || "AS" in params.analysis_method.split(",")) && params.pseudo_aligner in ["salmon", "kallisto"] ? "IsoformSwitchAnalyzeR (<a href=\"https://doi.org/10.1093/bioinformatics/btz247\" target=\"_blank\">Vitting-Seerup &amp; Sandelin, 2019</a>)," : "",
        (("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null)) || ("DEG" in params.analysis_method.split(",") && (params.msigdb_categories ? params.msigdb_categories.split(",").any { it in ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "H"] } : null)) ? "clusterProfiler (<a href=\"https://doi.org/10.1016/j.xinn.2024.100722\" target=\"_blank\">Yu, 2024</a>)," : "",
        "DEG" in params.analysis_method.split(",") && (params.enrichment_method ? "GSEA" in params.enrichment_method.split(",") : null) && (params.msigdb_categories ? params.msigdb_categories.split(",").any { it in ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "H"] } : null) ? "enrichplot (<a href=\"https://bioconductor.org/packages/enrichplot\" target=\"_blank\">Yu, 2025</a>)," : "",
        "DEG" in params.analysis_method.split(",") && (params.enrichment_method ? "GSEA" in params.enrichment_method.split(",") : null) && (params.msigdb_categories ? params.msigdb_categories.split(",").any { it in ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "H"] } : null) ? "msigdbr (<a href=\"https://doi.org/10.32614/cran.package.msigdbr\" target=\"_blank\">Dolgalev, 2025</a>)," : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) ? "purrr (<a href=\"https://doi.org/10.32614/cran.package.purrr\" target=\"_blank\">Wickham &amp; Henry, 2025</a>)," : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "human" ? "org.Hs.eg.db (<a href=\"https://doi.org/10.18129/B9.BIOC.ORG.HS.EG.DB\" target=\"_blank\">Carlson, 2025</a>)," : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "mouse" ? "org.Mm.eg.db (<a href=\"https://doi.org/10.18129/B9.BIOC.ORG.MM.EG.DB\" target=\"_blank\">Carlson, 2025</a>)," : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "rat" ? "org.Rn.eg.db (<a href=\"https://doi.org/10.18129/B9.BIOC.ORG.RN.EG.DB\" target=\"_blank\">Carlson, 2025</a>)," : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "yeast" ? "org.Sc.sgd.db (<a href=\"https://doi.org/10.18129/B9.BIOC.ORG.SC.SGD.DB\" target=\"_blank\">Carlson, 2025</a>)," : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "fruitfly" ? "org.Dm.eg.db (<a href=\"https://doi.org/10.18129/B9.BIOC.ORG.DM.EG.DB\" target=\"_blank\">Carlson, 2025</a>)," : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "zebrafish" ? "org.Dr.eg.db (<a href=\"https://doi.org/10.18129/B9.BIOC.ORG.DR.EG.DB\" target=\"_blank\">Carlson, 2025</a>)," : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "worm" ? "org.Ce.eg.db (<a href=\"https://doi.org/10.18129/B9.BIOC.ORG.CE.EG.DB\" target=\"_blank\">Carlson, 2025</a>)," : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "arabidopsis" ? "org.At.tair.db (<a href=\"https://doi.org/10.18129/B9.BIOC.ORG.AT.TAIR.DB\" target=\"_blank\">Carlson, 2025</a>)," : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "chicken" ? "org.Gg.eg.db (<a href=\"https://doi.org/10.18129/B9.BIOC.ORG.GG.EG.DB\" target=\"_blank\">Carlson, 2025</a>)," : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "cow" ? "org.Bt.eg.db (<a href=\"https://doi.org/10.18129/B9.BIOC.ORG.BT.EG.DB\" target=\"_blank\">Carlson, 2025</a>)," : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "pig" ? "org.Ss.eg.db (<a href=\"https://doi.org/10.18129/B9.BIOC.ORG.SS.EG.DB\" target=\"_blank\">Carlson, 2025</a>)," : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "dog" ? "org.Cf.eg.db (<a href=\"https://doi.org/10.18129/B9.BIOC.ORG.CF.EG.DB\" target=\"_blank\">Carlson, 2025</a>)," : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "monkey" ? "org.Mmu.eg.db (<a href=\"https://doi.org/10.18129/B9.BIOC.ORG.MMU.EG.DB\" target=\"_blank\">Carlson, 2025</a>)," : "",
    ].findAll { it }.join(' ').trim()

    return citation_text.replaceAll(", \\.", ".").replaceAll("\\. \\.", ".")
}

def toolBibliographyText() {
    def reference_text = [
        !params.skip_fastqc ? "<li>Andrews, S. (2010). FastQC: a quality control tool for high throughput sequence data. Retrieved from <a href=\"https://www.bioinformatics.babraham.ac.uk/projects/fastqc/\" target=\"_blank\">https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</a></li>" : "",
        !params.skip_multiqc ? "<li>Ewels, P., Magnusson, M., Lundin, S., &amp; Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. <em>Bioinformatics</em>, <em>32</em>(19), 3047–3048. doi: <a href=\"https://doi.org/10.1093/bioinformatics/btw354\" target=\"_blank\">10.1093/bioinformatics/btw354</a></li>" : "",
        "<li>Shen, W., Le, S., Li, Y., &amp; Hu, F. (2016). SeqKit: A cross-platform and ultrafast toolkit for FASTA/Q file manipulation. <em>PloS One</em>, <em>11</em>(10), e0163962. doi: <a href=\"https://doi.org/10.1371/journal.pone.0163962\" target=\"_blank\">10.1371/journal.pone.0163962</a></li>",
        params.aligner == "hisat2" ? "<li>Kim, D., Paggi, J. M., Park, C., Bennett, C., &amp; Salzberg, S. L. (2019). Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. <em>Nature Biotechnology</em>, <em>37</em>(8), 907–915. doi: <a href=\"https://doi.org/10.1038/s41587-019-0201-4\" target=\"_blank\">10.1038/s41587-019-0201-4</a></li>" : "",
        params.pseudo_aligner == "salmon" ? "<li>Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., &amp; Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. <em>Nature Methods</em>, <em>14</em>(4), 417–419. doi: <a href=\"https://doi.org/10.1038/nmeth.4197\" target=\"_blank\">10.1038/nmeth.4197</a></li>" : "",
        params.pseudo_aligner == "kallisto" ? "<li>Bray, N. L., Pimentel, H., Melsted, P., &amp; Pachter, L. (2016). Near-optimal probabilistic RNA-seq quantification. <em>Nature Biotechnology</em>, <em>34</em>(5), 525–527. doi: <a href=\"https://doi.org/10.1038/nbt.3519\" target=\"_blank\">10.1038/nbt.3519</a></li>" : "",
        params.with_umi ? "<li>Smith, T., Heger, A., &amp; Sudbery, I. (2017). UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. <em>Genome Research</em>, <em>27</em>(3), 491–499. doi: <a href=\"https://doi.org/10.1101/gr.209601.116\" target=\"_blank\">10.1101/gr.209601.116</a></li>" : "",
        !params.skip_trimming ? "<li>Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. <em>EMBnet.Journal</em>, <em>17</em>(1), 10. doi: <a href=\"https://doi.org/10.14806/ej.17.1.200\" target=\"_blank\">10.14806/ej.17.1.200</a></li>" : "",
        params.aligner == "hisat2" ? "<li>Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., … 1000 Genome Project Data Processing Subgroup. (2009). The Sequence Alignment/Map format and SAMtools. <em>Bioinformatics (Oxford, England)</em>, <em>25</em>(16), 2078–2079. doi: <a href=\"https://doi.org/10.1093/bioinformatics/btp352\" target=\"_blank\">10.1093/bioinformatics/btp352</a></li>" : "",
        params.aligner == "hisat2" && "DEG" in params.analysis_method.split(",") ? "<li>Liao, Y., Smyth, G. K., &amp; Shi, W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. <em>Bioinformatics (Oxford, England)</em>, <em>30</em>(7), 923–930. doi: <a href=\"https://doi.org/10.1093/bioinformatics/btt656\" target=\"_blank\">10.1093/bioinformatics/btt656</a></li>" : "",
        "AS" in params.analysis_method.split(",") && params.aligner == "hisat2" ? "<li>Wang, Y., Xie, Z., Kutschera, E., Adams, J. I., Kadash-Edmondson, K. E., &amp; Xing, Y. (2024). rMATS-turbo: an efficient and flexible computational tool for alternative splicing analysis of large-scale RNA-seq data. <em>Nature Protocols</em>, <em>19</em>(4), 1083–1104. doi: <a href=\"https://doi.org/10.1038/s41596-023-00944-2\" target=\"_blank\">10.1038/s41596-023-00944-2</a></li>" : "",
        "AS" in params.analysis_method.split(",") && params.aligner == "hisat2" ? "<li>Shieh, Kutschera, E., Tseng, Y.-T., DM_, Lin, I.-H., &amp; Samani, E. (2023). <em>Xinglab/rmats2sashimiplot: v3.0.0</em>. doi: <a href=\"https://doi.org/10.5281/ZENODO.10008656\" target=\"_blank\">10.5281/ZENODO.10008656</a></li>" : "",
        params.aligner == "hisat2" && (params.rseqc_modules ? params.rseqc_modules.split(",").any { it in ["bam_stat", "genebody_coverage", "infer_experiment", "inner_distance", "junction_annotation", "read_distribution", "read_duplication", "tin"] } : null) ? "<li>Wang, L., Wang, S., &amp; Li, W. (2012). RSeQC: quality control of RNA-seq experiments. <em>Bioinformatics (Oxford, England)</em>, <em>28</em>(16), 2184–2185. doi: <a href=\"https://doi.org/10.1093/bioinformatics/bts356\" target=\"_blank\">10.1093/bioinformatics/bts356</a></li>" : "",
        "DEU" in params.analysis_method.split(",") ? "<li>Anders, S., Reyes, A., &amp; Huber, W. (2012). Detecting differential usage of exons from RNA-seq data. <em>Genome Research</em>, <em>22</em>(10), 2008–2017. doi: <a href=\"https://doi.org/10.1101/gr.133744.111\" target=\"_blank\">10.1101/gr.133744.111</a></li>" : "",
        "DEG" in params.analysis_method.split(",") && (params.diffexpr_method == "limma" || params.diffexpr_method == "masigpro") ? "<li>Chen, Y., Chen, L., Lun, A. T. L., Baldoni, P. L., &amp; Smyth, G. K. (2025). edgeR v4: powerful differential analysis of sequencing data with expanded functionality and improved support for small counts and larger datasets. <em>Nucleic Acids Research</em>, <em>53</em>(2). doi: <a href=\"https://doi.org/10.1093/nar/gkaf018\" target=\"_blank\">10.1093/nar/gkaf018</a></li>" : "",
        "DEG" in params.analysis_method.split(",") && params.diffexpr_method == "limma" ? "<li>Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., &amp; Smyth, G. K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. <em>Nucleic Acids Research</em>, <em>43</em>(7), e47. doi: <a href=\"https://doi.org/10.1093/nar/gkv007\" target=\"_blank\">10.1093/nar/gkv007</a></li>" : "",
        ("DEG" in params.analysis_method.split(",") && params.diffexpr_method == "deseq2") || "WGCNA" in params.analysis_method.split(",") ? "<li>Love, M. I., Huber, W., &amp; Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. <em>Genome Biology</em>, <em>15</em>(12), 550. doi: <a href=\"https://doi.org/10.1186/s13059-014-0550-8\" target=\"_blank\">10.1186/s13059-014-0550-8</a></li>" : "",
        "DEG" in params.analysis_method.split(",") && params.diffexpr_method == "masigpro" ? "<li>Wickham, H. (2007). Reshaping Data with the reshape Package. <em>Journal of Statistical Software</em>, <em>21</em>(12). doi: <a href=\"https://doi.org/10.18637/jss.v021.i12\" target=\"_blank\">10.18637/jss.v021.i12</a></li>" : "",
        "DEG" in params.analysis_method.split(",") && params.diffexpr_method == "masigpro" ? "<li>Pedersen, T. L. (2019). patchwork: The Composer of Plots [Data set]. <em>CRAN: Contributed Packages</em>. doi: <a href=\"https://doi.org/10.32614/cran.package.patchwork\" target=\"_blank\">10.32614/cran.package.patchwork</a></li>" : "",
        "DEG" in params.analysis_method.split(",") && params.diffexpr_method == "masigpro" ? "<li>Scrucca, L., Fraley, C., Murphy, T. B., &amp; Raftery, A. E. (2023). Model-based clustering, classification, and density estimation using mclust in R. doi: <a href=\"https://doi.org/10.1201/9781003277965\" target=\"_blank\">10.1201/9781003277965</a></li>" : "",
        "DEG" in params.analysis_method.split(",") && params.pseudo_aligner in ["kallisto", "salmon"] ? "<li>Soneson, C., Love, M. I., &amp; Robinson, M. D. (2015). Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. <em>F1000Research</em>, <em>4</em>, 1521. doi: <a href=\"https://doi.org/10.12688/f1000research.7563.2\" target=\"_blank\">10.12688/f1000research.7563.2</a></li>" : "",
        "DEG" in params.analysis_method.split(",") && params.pseudo_aligner in ["kallisto", "salmon"] ? "<li>Fischer, B., Smith, M., &amp; Pau, G. (2025). <em>rhdf5: R Interface to HDF5</em>. Retrieved from <a href=\"https://bioconductor.org/packages/rhdf5\" target=\"_blank\">https://bioconductor.org/packages/rhdf5</a></li>" : "",
        "DEG" in params.analysis_method.split(",") && params.pseudo_aligner in ["kallisto", "salmon"] ? "<li>Ooms, J. (2014). The jsonlite package: A practical and consistent mapping between JSON data and R objects. doi: <a href=\"https://doi.org/10.48550/ARXIV.1403.2805\" target=\"_blank\">10.48550/ARXIV.1403.2805</a></li>" : "",
        "DEG" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || (("DIU" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null)) ? "<li>Durinck, S., Spellman, P. T., Birney, E., &amp; Huber, W. (2009). Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt. <em>Nature Protocols</em>, <em>4</em>(8), 1184–1191. doi: <a href=\"https://doi.org/10.1038/nprot.2009.97\" target=\"_blank\">10.1038/nprot.2009.97</a></li>" : "",
        "DEG" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") ? "<li>Kolde, R. (2025). <em>pheatmap: Pretty Heatmaps</em>. Retrieved from <a href=\"https://github.com/raivokolde/pheatmap\" target=\"_blank\">https://github.com/raivokolde/pheatmap</a></li>" : "",
        "DEG" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") ? "<li>Blighe, K., Rana, S., &amp; Lewis, M. (2019). <em>EnhancedVolcano: Publication-ready volcano plots with enhanced colouring and labeling</em> (R package version 1.0). Retrieved from <a href=\"https://github.com/kevinblighe/EnhancedVolcano\" target=\"_blank\">https://github.com/kevinblighe/EnhancedVolcano</a></li>" : "",
        "DEG" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") ? "<li>Neuwirth, E. (2002). RColorBrewer: ColorBrewer Palettes [Data set]. <em>CRAN: Contributed Packages</em>. doi: <a href=\"https://doi.org/10.32614/cran.package.rcolorbrewer\" target=\"_blank\">10.32614/cran.package.rcolorbrewer</a></li>" : "",
        "DEG" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || ("AS" in params.analysis_method.split(",") && params.pseudo_aligner in ["salmon", "kallisto"]) ? "<li>Wickham, H. (2016). <em>ggplot2: Elegant Graphics for Data Analysis</em>. Retrieved from <a href=\"https://ggplot2.tidyverse.org\" target=\"_blank\">https://ggplot2.tidyverse.org</a></li>" : "",
        "DEG" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || (("DIU" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null)) ? "<li>Müller, K., &amp; Wickham, H. (2016). tibble: Simple Data Frames [Data set]. <em>CRAN: Contributed Packages</em>. doi: <a href=\"https://doi.org/10.32614/cran.package.tibble\" target=\"_blank\">10.32614/cran.package.tibble</a></li>" : "",
        "DEG" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") ? "<li>Wickham, H., Vaughan, D., &amp; Girlich, M. (2014). tidyr: Tidy Messy Data [Data set]. <em>CRAN: Contributed Packages</em>. doi: <a href=\"https://doi.org/10.32614/cran.package.tidyr\" target=\"_blank\">10.32614/cran.package.tidyr</a></li>" : "",
        "DEG" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || (("DIU" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null)) ? "<li>Wickham, H., François, R., Henry, L., Müller, K., &amp; Vaughan, D. (2014). dplyr: A Grammar of Data Manipulation [Data set]. <em>CRAN: Contributed Packages</em>. doi: <a href=\"https://doi.org/10.32614/cran.package.dplyr\" target=\"_blank\">10.32614/cran.package.dplyr</a></li>" : "",
        "DEG" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") ? "<li>Wilke, C. (2025). cowplot: Streamlined Plot Theme and Plot Annotations for 'ggplot2' [Data set]. <em>CRAN: Contributed Packages</em>. doi: <a href=\"https://doi.org/10.32614/cran.package.cowplot\" target=\"_blank\">10.32614/cran.package.cowplot</a></li>" : "",
        "WGCNA" in params.analysis_method.split(",") ? "<li>Langfelder, P., &amp; Horvath, S. (2008). WGCNA: an R package for weighted correlation network analysis. <em>BMC Bioinformatics</em>, <em>9</em>(1), 559. doi: <a href=\"https://doi.org/10.1186/1471-2105-9-559\" target=\"_blank\">10.1186/1471-2105-9-559</a></li>" : "",
        "WGCNA" in params.analysis_method.split(",") ? "<li>Szklarczyk, D., Kirsch, R., Koutrouli, M., Nastou, K., Mehryary, F., Hachilif, R., … von Mering, C. (2023). The STRING database in 2023: protein-protein association networks and functional enrichment analyses for any sequenced genome of interest. <em>Nucleic Acids Research</em>, <em>51</em>(D1), D638–D646. doi: <a href=\"https://doi.org/10.1093/nar/gkac1000\" target=\"_blank\">10.1093/nar/gkac1000</a></li>" : "",
        ("DIU" in params.analysis_method.split(",") || "AS" in params.analysis_method.split(",")) && params.pseudo_aligner in ["salmon", "kallisto"] ? "<li>Vitting-Seerup, K., &amp; Sandelin, A. (2019). IsoformSwitchAnalyzeR: analysis of changes in genome-wide patterns of alternative splicing and its functional consequences. <em>Bioinformatics (Oxford, England)</em>, <em>35</em>(21), 4469–4471. doi: <a href=\"https://doi.org/10.1093/bioinformatics/btz247\" target=\"_blank\">10.1093/bioinformatics/btz247</a></li>" : "",
        (("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null)) || ("DEG" in params.analysis_method.split(",") && (params.msigdb_categories ? params.msigdb_categories.split(",").any { it in ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "H"] } : null)) ? "<li>Yu, G. (2024). Thirteen years of clusterProfiler. <em>Innovation (Cambridge (Mass.))</em>, <em>5</em>(6), 100722. doi: <a href=\"https://doi.org/10.1016/j.xinn.2024.100722\" target=\"_blank\">10.1016/j.xinn.2024.100722</a></li>" : "",
        "DEG" in params.analysis_method.split(",") && (params.enrichment_method ? "GSEA" in params.enrichment_method.split(",") : null) && (params.msigdb_categories ? params.msigdb_categories.split(",").any { it in ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "H"] } : null) ? "<li>Yu, G. (2025). <em>enrichplot: Visualization of Functional Enrichment Result</em> (R package version 1.30.4). Retrieved from <a href=\"https://bioconductor.org/packages/enrichplot\" target=\"_blank\">https://bioconductor.org/packages/enrichplot</a></li>" : "",
        "DEG" in params.analysis_method.split(",") && (params.enrichment_method ? "GSEA" in params.enrichment_method.split(",") : null) && (params.msigdb_categories ? params.msigdb_categories.split(",").any { it in ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "H"] } : null) ? "<li>Dolgalev, I. (2025). Msigdbr: MSigDB gene sets for multiple organisms in a tidy data format [Data set]. <em>CRAN: Contributed Packages</em>. doi: <a href=\"https://doi.org/10.32614/cran.package.msigdbr\" target=\"_blank\">10.32614/cran.package.msigdbr</a></li>" : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) ? "<li>Wickham, H., &amp; Henry, L. (2025). purrr: Functional Programming Tools [Data set]. <em>CRAN: Contributed Packages</em>. doi: <a href=\"https://doi.org/10.32614/cran.package.purrr\" target=\"_blank\">10.32614/cran.package.purrr</a></li>" : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "human" ? "<li>Carlson, M. (2025). <em>org.Hs.eg.db</em>. doi: <a href=\"https://doi.org/10.18129/B9.BIOC.ORG.HS.EG.DB\" target=\"_blank\">10.18129/B9.BIOC.ORG.HS.EG.DB</a></li>" : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "mouse" ? "<li>Carlson, M. (2025). <em>org.Mm.eg.db</em>. doi: <a href=\"https://doi.org/10.18129/B9.BIOC.ORG.MM.EG.DB\" target=\"_blank\">10.18129/B9.BIOC.ORG.MM.EG.DB</a></li>" : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "rat" ? "<li>Carlson, M. (2025). <em>org.Rn.eg.db</em>. doi: <a href=\"https://doi.org/10.18129/B9.BIOC.ORG.RN.EG.DB\" target=\"_blank\">10.18129/B9.BIOC.ORG.RN.EG.DB</a></li>" : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "yeast" ? "<li>Carlson, M. (2025). <em>org.Sc.sgd.db</em>. doi: <a href=\"https://doi.org/10.18129/B9.BIOC.ORG.SC.SGD.DB\" target=\"_blank\">10.18129/B9.BIOC.ORG.SC.SGD.DB</a></li>" : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "fruitfly" ? "<li>Carlson, M. (2025). <em>org.Dm.eg.db</em>. doi: <a href=\"https://doi.org/10.18129/B9.BIOC.ORG.DM.EG.DB\" target=\"_blank\">10.18129/B9.BIOC.ORG.DM.EG.DB</a></li>" : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "zebrafish" ? "<li>Carlson, M. (2025). <em>org.Dr.eg.db</em>. doi: <a href=\"https://doi.org/10.18129/B9.BIOC.ORG.DR.EG.DB\" target=\"_blank\">10.18129/B9.BIOC.ORG.DR.EG.DB</a></li>" : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "worm" ? "<li>Carlson, M. (2025). <em>org.Ce.eg.db</em>. doi: <a href=\"https://doi.org/10.18129/B9.BIOC.ORG.CE.EG.DB\" target=\"_blank\">10.18129/B9.BIOC.ORG.CE.EG.DB</a></li>" : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "arabidopsis" ? "<li>Carlson, M. (2025). <em>org.at.tair.db</em>. doi: <a href=\"https://doi.org/10.18129/B9.BIOC.ORG.AT.TAIR.DB\" target=\"_blank\">10.18129/B9.BIOC.ORG.AT.TAIR.DB</a></li>" : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "chicken" ? "<li>Carlson, M. (2025). <em>org.Gg.eg.db</em>. doi: <a href=\"https://doi.org/10.18129/B9.BIOC.ORG.GG.EG.DB\" target=\"_blank\">10.18129/B9.BIOC.ORG.GG.EG.DB</a></li>" : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "cow" ? "<li>Carlson, M. (2025). <em>org.Bt.eg.db</em>. doi: <a href=\"https://doi.org/10.18129/B9.BIOC.ORG.BT.EG.DB\" target=\"_blank\">10.18129/B9.BIOC.ORG.BT.EG.DB</a></li>" : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "pig" ? "<li>Carlson, M. (2025). <em>org.Ss.eg.db</em>. doi: <a href=\"https://doi.org/10.18129/B9.BIOC.ORG.SS.EG.DB\" target=\"_blank\">10.18129/B9.BIOC.ORG.SS.EG.DB</a></li>" : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "dog" ? "<li>Carlson, M. (2025). <em>org.Cf.eg.db</em>. doi: <a href=\"https://doi.org/10.18129/B9.BIOC.ORG.CF.EG.DB\" target=\"_blank\">10.18129/B9.BIOC.ORG.CF.EG.DB</a></li>" : "",
        ("DEG" in params.analysis_method.split(",") || "DIU" in params.analysis_method.split(",") || "WGCNA" in params.analysis_method.split(",") || "DEU" in params.analysis_method.split(",")) && (params.enrichment_method ? params.enrichment_method.split(",").any { it in ["KEGG", "GO"] } : null) && params.species == "monkey" ? "<li>Carlson, M. (2025). <em>org.Mmu.eg.db</em>. doi: <a href=\"https://doi.org/10.18129/B9.BIOC.ORG.MMU.EG.DB\" target=\"_blank\">10.18129/B9.BIOC.ORG.MMU.EG.DB</a></li>" : "",
    ].findAll { it }.join(' ')

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert to a named map so can be used as with familiar NXF ${workflow} variable syntax
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    }
    else {
        meta["doi_text"] = ""
    }
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used.</li>"

    // Tool references
    meta["tool_citations"] = toolCitationText()
    meta["tool_bibliography"] = toolBibliographyText()

    def methods_text = mqc_methods_yaml.text

    def engine = new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
