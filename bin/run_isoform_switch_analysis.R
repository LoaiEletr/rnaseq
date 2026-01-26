#!/usr/bin/env Rscript

## =============================================================================
## Script:  run_isoform_switch_analysis.R
## Purpose: Isoform switching and alternative splicing analysis using IsoformSwitchAnalyzeR
## Usage:   Rscript run_isoform_switch_analysis.R <quant_dir1> <quant_dir2> ... <samplesheet>
##          <quant_type> <gtf_file> <transcript_fasta> <method> <dIF_cutoff> <pval_threshold> <top_n>
## Example: Rscript run_isoform_switch_analysis.R quant1 quant2 samples.csv salmon ref.gtf
##          transcripts.fa DIU,AS 0.1 0.05 10
## =============================================================================

## Load required packages ------------------------------------------------------
library(IsoformSwitchAnalyzeR)
library(ggplot2)
library(MASS)
library(rhdf5)

## =============================================================================
## Parse command line arguments
## =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 9) {
    stop("Usage: run_isoform_switch_analysis.R <quant_dir1> <quant_dir2> ... <samplesheet> ",
         "<quant_type> <gtf_file> <transcript_fasta> <method> <dIF_cutoff> <pval_threshold> <top_n>",
         call. = FALSE)
}

## Calculate number of quantification directories ------------------------------
nQuantDirs <- length(args) - 8

## Extract parameters with descriptive names -----------------------------------
quantDirs <- args[1:nQuantDirs]
samplesheetFile <- args[nQuantDirs + 1]
quantType <- args[nQuantDirs + 2]
gtfFile <- args[nQuantDirs + 3]
transcriptFasta <- args[nQuantDirs + 4]
methodInput <- args[nQuantDirs + 5]
dIfCutoff <- as.numeric(args[nQuantDirs + 6])
pvalueThreshold <- as.numeric(args[nQuantDirs + 7])
topNIsoforms <- as.numeric(args[nQuantDirs + 8])

## =============================================================================
## Validate quantification type
## =============================================================================

## Set quantification file paths -----------------------------------------------
if (quantType == "salmon") {
    quantFiles <- file.path(quantDirs, "quant.sf")
} else if (quantType == "kallisto") {
    quantFiles <- file.path(quantDirs, "abundance.tsv")
} else {
    stop("Error: quant_type must be either 'salmon' or 'kallisto'", call. = FALSE)
}

## Convert method input to vector ----------------------------------------------
analysisMethods <- unlist(strsplit(methodInput, ","))

## =============================================================================
## Load and prepare sample information
## =============================================================================

## Load sample information -----------------------------------------------------
samplesheet <- read.csv(samplesheetFile, stringsAsFactors = FALSE)

## Subset samplesheet to required columns --------------------------------------
samplesheet <- samplesheet[, c("sample_id", "condition"), drop = FALSE]
colnames(samplesheet)[colnames(samplesheet) == "sample_id"] <- "sampleID"

## =============================================================================
## Import isoform expression data
## =============================================================================

cat("Importing isoform expression data...\n")
isoformData <- importIsoformExpression(sampleVector = quantFiles)

## Rename columns to match sample IDs -----------------------------------------
colnames(isoformData$abundance) <- c("isoform_id", basename(quantDirs))
colnames(isoformData$counts) <- c("isoform_id", basename(quantDirs))

## Reorder columns to match samplesheet ---------------------------------------
isoformData$counts <- isoformData$counts[, c("isoform_id", samplesheet$sampleID)]
isoformData$abundance <- isoformData$abundance[, c("isoform_id", samplesheet$sampleID)]

## =============================================================================
## Create output directory
## =============================================================================

outputDir <- "isoform_output"
dir.create(outputDir, showWarnings = FALSE)

## =============================================================================
## Import data into IsoformSwitchAnalyzeR object with error handling
## =============================================================================

cat("Creating IsoformSwitchAnalyzeR object...\n")

tryCatch({
    isoformSwitchList <- importRdata(
        isoformCountMatrix = isoformData$counts,
        isoformRepExpression = isoformData$abundance,
        designMatrix = samplesheet,
        removeNonConvensionalChr = TRUE,
        addAnnotatedORFs = TRUE,
        isoformExonAnnoation = gtfFile,
        isoformNtFasta = transcriptFasta,
        showProgress = TRUE
    )
}, error = function(error) {

    ## Save detailed error message ---------------------------------------------
    errorText <- paste(
        "Analysis stopped: importRdata failed",
        paste("Error:", error$message),
        paste("Time:", Sys.time()),
        "",
        "Common reasons for this error:",
        "1. GTF and quantification IDs don't match",
        "2. Too few overlapping transcripts",
        "3. GTF file format issues",
        "4. Non-matching annotation files",
        sep = "\n"
    )

    writeLines(errorText, file.path(outputDir, "import_error.txt"))

    message("Error: ", error$message)
    quit(save = "no", status = 0)
})

## =============================================================================
## Helper function to save tables
## =============================================================================

saveTableCsv <- function(dataTable, fileName, outputPath) {
    csvPath <- file.path(outputPath, "csv", fileName)
    write.csv(dataTable, csvPath, row.names = FALSE)
}

## =============================================================================
## Differential Isoform Usage (DIU) analysis
## =============================================================================

if ("DIU" %in% analysisMethods) {
    cat("\nPerforming Differential Isoform Usage (DIU) analysis...\n")

    tryCatch({
        switchAnalysisList <- isoformSwitchAnalysisCombined(
            switchAnalyzeRlist = isoformSwitchList,
            pathToOutput = outputDir
        )
    }, error = function(error) {
        message("DIU analysis stopped due to: ", error$message)
        quit(status = 0)
    })

    ## DIU Volcano plot --------------------------------------------------------
    diuDir <- file.path(outputDir, "DIU")
    dir.create(diuDir, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(diuDir, "csv"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(diuDir, "top_switches"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(diuDir, "plots"), recursive = TRUE, showWarnings = FALSE)

    tryCatch({
        diuVolcanoPlot <- ggplot(
            data = switchAnalysisList$isoformFeatures,
            aes(x = dIF, y = -log10(isoform_switch_q_value))
        ) +
            geom_point(
                aes(color = abs(dIF) > dIfCutoff & isoform_switch_q_value < pvalueThreshold),
                size = 1
            ) +
            geom_hline(yintercept = -log10(pvalueThreshold), linetype = 'dashed') +
            geom_vline(xintercept = c(-dIfCutoff, dIfCutoff), linetype = 'dashed') +
            scale_color_manual('Significant\nIsoform Switch', values = c('black', 'red')) +
            labs(x = 'dIF', y = '-Log10 (Isoform Switch Q Value)') +
            theme_bw()

        ggsave(
            file.path(diuDir, "plots/DIU_volcano_plot.pdf"),
            diuVolcanoPlot,
            width = 10,
            height = 8
        )
    }, error = function(error) {
        message("DIU volcano plot generation failed: ", error$message)
    })

    ## Identify significant DIU genes ------------------------------------------
    isoformFeatures <- switchAnalysisList$isoformFeatures
    significantDiuGenes <- isoformFeatures[
        abs(isoformFeatures$dIF) > dIfCutoff &
        isoformFeatures$isoform_switch_q_value < pvalueThreshold,
    ]

    saveTableCsv(significantDiuGenes, "significant_DIU_genes.csv", diuDir)

    significantDiuGeneIds <- unique(
        significantDiuGenes[
            significantDiuGenes$switchConsequencesGene == TRUE,
            "gene_id"
        ]
    )

    saveRDS(
        significantDiuGeneIds,
        file.path(outputDir, "significant_DIU_gene_ids.rds")
    )

    ## Extract top switches ----------------------------------------------------
    tryCatch({
        topSwitchesTable <- extractTopSwitches(
            switchAnalysisList,
            filterForConsequences = TRUE,
            alpha = pvalueThreshold,
            dIFcutoff = dIfCutoff,
            n = topNIsoforms,
        )

        saveTableCsv(
            topSwitchesTable,
            paste0("top_", topNIsoforms, "_switches.csv"),
            diuDir
        )
    }, error = function(error) {
        message("Failed to extract top switches: ", error$message)
    })

    ## Generate plots for top switches -----------------------------------------
    tryCatch({
        switchPlotTopSwitches(
            switchAnalysisList,
            alpha = pvalueThreshold,
            dIFcutoff = dIfCutoff,
            n = topNIsoforms,
            pathToOutput = file.path(diuDir, "top_switches"),
        )
    }, error = function(error) {
        message("Failed to generate top switch plots: ", error$message)
    })

    ## Consequence analysis ----------------------------------------------------
    cat("Adding consequence analysis to DIU output...\n")

    ### Consequence enrichment plot
    tryCatch({
        consequenceEnrichmentPlot <- extractConsequenceEnrichment(
            switchAnalysisList,
            alpha = pvalueThreshold,
            dIFcutoff = dIfCutoff,
            consequencesToAnalyze = 'all',
            analysisOppositeConsequence = TRUE,
            localTheme = theme_bw(base_size = 14),
            returnResult = FALSE
        )
        ggsave(
            file.path(diuDir, "plots/consequence_enrichment_plot.pdf"),
            consequenceEnrichmentPlot,
            height = 9,
            width = 10
        )
    }, error = function(error) {
        message("Failed to create consequence enrichment plot: ", error$message)
    })

    ### Consequence enrichment data
    tryCatch({
        consequenceEnrichmentData <- extractConsequenceEnrichment(
            switchAnalysisList,
            alpha = pvalueThreshold,
            dIFcutoff = dIfCutoff,
            consequencesToAnalyze = 'all',
            analysisOppositeConsequence = TRUE,
            returnResult = TRUE
        )
        saveTableCsv(consequenceEnrichmentData, "consequence_enrichment_data.csv", diuDir)
    }, error = function(error) {
        message("Failed to extract consequence enrichment data: ", error$message)
    })

    ### Consequence genomewide plot
    tryCatch({
        consequenceGenomewidePlot <- extractConsequenceGenomeWide(
            switchAnalysisList,
            alpha = pvalueThreshold,
            dIFcutoff = dIfCutoff,
            returnResult = FALSE
        )
        ggsave(
            file.path(diuDir, "plots/consequence_genomewide_plot.pdf"),
            consequenceGenomewidePlot,
            height = 9,
            width = 20
        )
    }, error = function(error) {
        message("Failed to create consequence genomewide plot: ", error$message)
    })

    ### Consequence genomewide data
    tryCatch({
        consequenceGenomewideData <- extractConsequenceGenomeWide(
            switchAnalysisList,
            alpha = pvalueThreshold,
            dIFcutoff = dIfCutoff,
            returnResult = TRUE
        )
        saveTableCsv(consequenceGenomewideData, "consequence_genomewide_data.csv", diuDir)
    }, error = function(error) {
        message("Failed to extract consequence genomewide data: ", error$message)
    })

    ### Consequence summary plot
    tryCatch({
        consequenceSummaryPlot <- extractConsequenceSummary(
            switchAnalysisList,
            alpha = pvalueThreshold,
            dIFcutoff = dIfCutoff,
            returnResult = FALSE
        )
        ggsave(
            file.path(diuDir, "plots/consequence_summary_plot.pdf"),
            consequenceSummaryPlot,
            height = 9,
            width = 10
        )
    }, error = function(error) {
        message("Failed to create consequence summary plot: ", error$message)
    })

    ### Consequence summary data
    tryCatch({
        consequenceSummaryData <- extractConsequenceSummary(
            switchAnalysisList,
            alpha = pvalueThreshold,
            dIFcutoff = dIfCutoff,
            returnResult = TRUE
        )
        saveTableCsv(consequenceSummaryData, "consequence_summary_data.csv", diuDir)
    }, error = function(error) {
        message("Failed to extract consequence summary data: ", error$message)
    })

    ### Switch summary data
    tryCatch({
        switchSummaryData <- extractSwitchSummary(
        switchAnalysisList,
        alpha = pvalueThreshold,
        dIFcutoff = dIfCutoff,
        )
        saveTableCsv(switchSummaryData, "switch_summary_data.csv", asDir)
    }, error = function(error) {
        message("Failed to extract switch summary data: ", error$message)
    })
}

## =============================================================================
## Alternative Splicing (AS) analysis
## =============================================================================

if ("AS" %in% analysisMethods) {
    cat("\nPerforming Alternative Splicing (AS) analysis...\n")

    ## Create output directory for AS results ----------------------------------
    asDir <- file.path(outputDir, "AS")
    dir.create(asDir, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(asDir, "csv"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(asDir, "plots"), recursive = TRUE, showWarnings = FALSE)

    ## Use switchList from DIU analysis or re-run if needed --------------------
    if (!exists("switchAnalysisList")) {
        tryCatch({
            switchAnalysisList <- isoformSwitchAnalysisCombined(
                switchAnalyzeRlist = isoformSwitchList,
                pathToOutput = outputDir
            )
        }, error = function(error) {
            message("AS analysis stopped due to: ", error$message)
            quit(status = 0)
        })
    }

    ## Splicing enrichment analysis --------------------------------------------
    tryCatch({
        splicingEnrichmentPlot <- extractSplicingEnrichment(
            switchAnalysisList,
            alpha = pvalueThreshold,
            dIFcutoff = dIfCutoff,
            returnResult = FALSE
        )
        ggsave(
            file.path(asDir, "plots/splicing_enrichment_plot.pdf"),
            splicingEnrichmentPlot,
            height = 9,
            width = 10
        )
    }, error = function(error) {
        message("Failed to create splicing enrichment plot: ", error$message)
    })

    ### Splicing enrichment data
    tryCatch({
        splicingEnrichmentData <- extractSplicingEnrichment(
            switchAnalysisList,
            alpha = pvalueThreshold,
            dIFcutoff = dIfCutoff,
            returnResult = TRUE
        )
        saveTableCsv(splicingEnrichmentData, "splicing_enrichment_data.csv", asDir)
    }, error = function(error) {
        message("Failed to extract splicing enrichment data: ", error$message)
    })

    ## Splicing summary plot ---------------------------------------------------
    tryCatch({
        splicingSummaryPlot <- extractSplicingSummary(
            switchAnalysisList,
            alpha = pvalueThreshold,
            dIFcutoff = dIfCutoff,
            returnResult = FALSE
        )
        ggsave(
            file.path(asDir, "plots/splicing_summary_plot.pdf"),
            splicingSummaryPlot,
            height = 9,
            width = 10
        )
    }, error = function(error) {
        message("Failed to create splicing summary plot: ", error$message)
    })

    ### Splicing summary data
    tryCatch({
        splicingSummaryData <- extractSplicingSummary(
            switchAnalysisList,
            alpha = pvalueThreshold,
            dIFcutoff = dIfCutoff,
            returnResult = TRUE
        )
        saveTableCsv(splicingSummaryData, "splicing_summary_data.csv", asDir)
    }, error = function(error) {
        message("Failed to extract splicing summary data: ", error$message)
    })

    ## Splicing genomewide plot ------------------------------------------------
    tryCatch({
        splicingGenomewidePlot <- extractSplicingGenomeWide(
            switchAnalysisList,
            alpha = pvalueThreshold,
            dIFcutoff = dIfCutoff,
            returnResult = FALSE
        )
        ggsave(
            file.path(asDir, "plots/splicing_genomewide_plot.pdf"),
            splicingGenomewidePlot,
            height = 9,
            width = 20
        )
    }, error = function(error) {
        message("Failed to create splicing genomewide plot: ", error$message)
    })

    ### Splicing genomewide data
    tryCatch({
        splicingGenomewideData <- extractSplicingGenomeWide(
            switchAnalysisList,
            alpha = pvalueThreshold,
            dIFcutoff = dIfCutoff,
            returnResult = TRUE
        )
        saveTableCsv(splicingGenomewideData, "splicing_genomewide_data.csv", asDir)
    }, error = function(error) {
        message("Failed to extract splicing genomewide data: ", error$message)
    })

    ## Switch summary data -----------------------------------------------------
    tryCatch({
        switchSummaryData <- extractSwitchSummary(
            switchAnalysisList,
            alpha = pvalueThreshold,
            dIFcutoff = dIfCutoff,
        )
        saveTableCsv(switchSummaryData, "switch_summary_data.csv", asDir)
    }, error = function(error) {
        message("Failed to extract switch summary data: ", error$message)
    })

    cat("AS analysis completed. Results saved in:", asDir, "\n")
}

## =============================================================================
## Save final objects
## =============================================================================

saveRDS(isoformSwitchList, file.path(outputDir, "switch_analyzeR_list.rds"))
saveRDS(isoformData$abundance, file.path(outputDir, "abundance.rds"))
saveRDS(isoformData$counts, file.path(outputDir, "counts.rds"))

## Remove corrupted/unreadable PDF files generated during analysis
file.remove(file.path(outputDir, c("switch_consequences_enrichment.pdf", "splicing_enrichment.pdf")))
