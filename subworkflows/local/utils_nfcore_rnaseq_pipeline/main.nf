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
    aligner //  string: Aligner to use (hisat2/salmon/kallisto)
    pseudo_aligner //  string: Pseudo-aligner to use (salmon/kallisto)
    analysis_method //  string: Analysis methods to run (DEG,WGCNA,DIU,AS)
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
\033[0;35m  LoaiEletr/rnaseq ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/', '')}" }.join("\n")}${workflow.manifest.doi ? "\n" : ""}

* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/LoaiEletr/rnaseq/blob/master/CITATIONS.md
"""
    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --genome GRCh38 --outdir <OUTDIR>"

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
    ch_samplesheet = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            // Extract all columns from YOUR custom samplesheet
            def meta = row.sample_id
            def fastq_1 = row.fastq_1
            def fastq_2 = row.fastq_2
            def condition = row.condition
            def lib_type = row.lib_type
            def sequencer = row.sequencer

            // Create the complete meta map matching your workflow
            def meta_map = [
                id: meta,
                condition: condition,
                lib_type: lib_type,
                sequencer: sequencer,
            ]

            // Handle single-end vs paired-end
            def fastqs
            if (!fastq_2 || fastq_2.trim().isEmpty()) {
                meta_map.single_end = true
                fastqs = [fastq_1]
            }
            else {
                meta_map.single_end = false
                fastqs = [fastq_1, fastq_2]
            }

            // Return exactly what your workflow expects
            return [meta_map, fastqs]
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .map { meta, fastqs ->
            return [meta, fastqs.flatten()]
        }
        .set { ch_samplesheet }

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

    // Validate aligner selection
    if (!["hisat2"].contains(params.aligner) && !["salmon", "kallisto"].contains(params.pseudo_aligner)) {
        error(
            """
        ❌ INVALID ALIGNER SELECTION

        Valid aligners are:
          --aligner 'hisat2'            (for HISAT2 alignment)
          --pseudo_aligner 'salmon'     (for Salmon quantification)
          --pseudo_aligner 'kallisto'   (for Kallisto quantification)

        You specified: --aligner '${params.aligner}' --pseudo_aligner '${params.pseudo_aligner}'
        """
        )
    }

    // Validate analysis methods
    def valid_methods = ["DEG", "WGCNA", "DIU", "AS"]
    def requested_methods = params.analysis_method ? params.analysis_method.split(",").collect { it.trim() } : []
    def invalid_methods = requested_methods.findAll { !valid_methods.contains(it) }

    if (invalid_methods) {
        error(
            """
        ❌ INVALID ANALYSIS METHOD(S): ${invalid_methods.join(", ")}

        Valid analysis methods are:
          • DEG   - Differential Expression (limma/DESeq2/edgeR)
          • WGCNA - Co-expression Network Analysis
          • DIU   - Differential Isoform Usage
          • AS    - Alternative Splicing (rMATS)

        Example usage:
          --analysis_method DEG,WGCNA   # Run DE and co-expression
          --analysis_method DIU,AS      # Run isoform and splicing
        """
        )
    }

    // Validate combination of aligner and analysis methods
    if (params.aligner == "hisat2" && "DIU" in requested_methods) {
        log.warn("⚠️  DIU analysis requires pseudo-aligner (salmon/kallisto), but HISAT2 is selected. Consider using --pseudo_aligner instead.")
    }

    // Check for required files based on analysis methods
    if (["DEG", "WGCNA"].any { it in requested_methods } && !params.genome) {
        log.warn("⚠️  DEG/WGCNA analysis recommended with genome annotation. Consider providing --genome parameter.")
    }
}

//
// Validate channels from input samplesheet - UPDATED FOR YOUR FORMAT
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect { meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    // Additional validation for YOUR custom columns
    if (metas.size() > 1) {
        // Check condition for differential expression
        def conditions = metas.collect { it.condition }.unique()
        if (conditions.size() < 2 && "DEG" in params.analysis_method?.split(",")) {
            log.warn("⚠️  Differential expression analysis requires at least 2 conditions. Found: ${conditions}")
        }

        // Check lib_type consistency
        def lib_types = metas.collect { it.lib_type }.unique()
        if (lib_types.size() > 1) {
            log.warn("⚠️  Multiple library types detected: ${lib_types}. Ensure this is intentional.")
        }

        // Validate lib_type values
        def valid_lib_types = ["forward", "reverse", "auto"]
        def invalid_lib_types = lib_types.findAll { !valid_lib_types.contains(it) }
        if (invalid_lib_types) {
            error("Invalid lib_type values: ${invalid_lib_types}. Valid values are: forward, reverse, auto")
        }

        // Check sequencer information
        def sequencers = metas.collect { it.sequencer }.unique()
        if (sequencers.size() > 1) {
            log.info("Multiple sequencers detected: ${sequencers}. Batch effects may need consideration.")
        }
    }

    return [metas[0], fastqs]
}
//
// Get attribute from genome config file e.g. fasta, gtf
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[params.genome].containsKey(attribute)) {
            return params.genomes[params.genome][attribute]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" + "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" + "  Currently, the available genome keys are:\n" + "  ${params.genomes.keySet().join(", ")}\n" + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC - RNA-seq specific
//
def toolCitationText() {
    def citation_text = [
        "Tools used in the workflow included:",
        !params.skip_fastqc ? "FastQC (Andrews 2010)," : "",
        !params.skip_multiqc ? "MultiQC (Ewels et al. 2016)," : "",
        params.aligner == "hisat2" ? "HISAT2 (Kim et al. 2019)," : "",
        params.pseudo_aligner == "salmon" ? "Salmon (Patro et al. 2017)," : "",
        params.pseudo_aligner == "kallisto" ? "Kallisto (Bray et al. 2016)," : "",
        params.with_umi ? "UMI-tools (Smith et al. 2017)," : "",
        !params.skip_trimming ? "Cutadapt (Martin 2011)," : "",
        "SAMtools (Li et al. 2009),",
        "featureCounts (Liao et al. 2014),",
        "RSeQC (Wang et al. 2012),",
        "DEXSeq (Anders et al. 2012),",
        "rMATS (Shen et al. 2014),",
        "limma (Ritchie et al. 2015),",
        "DESeq2 (Love et al. 2014),",
        "edgeR (Robinson et al. 2010),",
        "WGCNA (Langfelder and Horvath 2008),",
        "IsoformSwitchAnalyzeR (Vitting-Seerup et al. 2019),",
        "SeqKit (Shen et al. 2016)",
    ].findAll { it }.join(' ').trim()

    return citation_text.replaceAll(", \\.", ".").replaceAll("\\. \\.", ".")
}

def toolBibliographyText() {
    def reference_text = [
        "<li>Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</li>",
        "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047-3048.</li>",
        params.aligner == "hisat2" ? "<li>Kim, D., Paggi, J. M., Park, C., Bennett, C., & Salzberg, S. L. (2019). Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nature Biotechnology, 37(8), 907-915.</li>" : "",
        params.pseudo_aligner == "salmon" ? "<li>Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods, 14(4), 417-419.</li>" : "",
        params.pseudo_aligner == "kallisto" ? "<li>Bray, N. L., Pimentel, H., Melsted, P., & Pachter, L. (2016). Near-optimal probabilistic RNA-seq quantification. Nature Biotechnology, 34(5), 525-527.</li>" : "",
        params.with_umi ? "<li>Smith, T., Heger, A., & Sudbery, I. (2017). UMI-tools: modelling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. PeerJ, 5, e8275.</li>" : "",
        !params.skip_trimming ? "<li>Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10-12.</li>" : "",
        "<li>Li, H., Handsaker, B., Wysoker, A., et al. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078-2079.</li>",
        "<li>Liao, Y., Smyth, G. K., & Shi, W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7), 923-930.</li>",
        "<li>Wang, L., Wang, S., & Li, W. (2012). RSeQC: quality control of RNA-seq experiments. Bioinformatics, 28(16), 2184-2185.</li>",
        "<li>Anders, S., Reyes, A., & Huber, W. (2012). Detecting differential usage of exons from RNA-seq data. Genome Research, 22(10), 2008-2017.</li>",
        "<li>Shen, S., Park, J. W., Lu, Z. X., et al. (2014). rMATS: robust and flexible detection of differential alternative splicing from replicate RNA-Seq data. Proceedings of the National Academy of Sciences, 111(51), E5593-E5601.</li>",
        "<li>Ritchie, M. E., Phipson, B., Wu, D., et al. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research, 43(7), e47.</li>",
        "<li>Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), 550.</li>",
        "<li>Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26(1), 139-140.</li>",
        "<li>Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 9, 559.</li>",
        "<li>Vitting-Seerup, K., Sandelin, A., & Waage, J. (2019). IsoformSwitchAnalyzeR: analysis of changes in genome-wide patterns of alternative splicing and its functional consequences. Bioinformatics, 35(21), 4469-4471.</li>",
        "<li>Shen, W., Le, S., Li, Y., & Hu, F. (2016). SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PLoS ONE, 11(10), e0163962.</li>",
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
