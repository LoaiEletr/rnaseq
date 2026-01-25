#!/usr/bin/env Rscript

## =============================================================================
## Script:  run_dexseq_deu.R
## Purpose: Differential exon usage analysis using DEXSeq
## Usage:   Rscript run_dexseq_deu.R <samplesheet_csv> <flattened_gff> <pval_threshold>
##          <lfc_threshold> <top_n> <min_exon_length>
## Example: Rscript run_dexseq_deu.R samples.csv flattened.gff 0.05 1.0 20 50
## =============================================================================

## Load required packages ------------------------------------------------------
library(DEXSeq)

## =============================================================================
## Parse command line arguments
## =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
    stop("Usage: run_dexseq_deu.R <samplesheet_csv> <flattened_gff> <pvalue_threshold> ",
         "<logfc_threshold> <top_n_genes> <min_exon_length>",
         call. = FALSE)
}

## Extract parameters with descriptive names -----------------------------------
samplesheetFile <- args[1]
flattenedGffFile <- args[2]
pvalueThreshold <- as.numeric(args[3])
logFcThreshold <- as.numeric(args[4])
topNGenes <- as.numeric(args[5])
minExonLength <- as.numeric(args[6])

## =============================================================================
## Load and validate input data
## =============================================================================

## Load sample information -----------------------------------------------------
samplesheet <- read.csv(samplesheetFile, stringsAsFactors = FALSE)

## Subset samplesheet to required columns --------------------------------------
samplesheet <- samplesheet[, c("sample_id", "condition"), drop = FALSE]

## Validate required columns ---------------------------------------------------
if (!"condition" %in% colnames(samplesheet)) {
    stop("Column 'condition' not found in samplesheet", call. = FALSE)
}

## Ensure condition is a factor ------------------------------------------------
samplesheet$condition <- factor(samplesheet$condition)
rownames(samplesheet) <- samplesheet$sample_id

## =============================================================================
## Helper functions
## =============================================================================

## Clean DEXSeq input files ----------------------------------------------------
cleanDexseqCountFile <- function(inputFilePath) {
    fileLines <- readLines(inputFilePath)

    ## Remove lines starting with '__'
    fileLines <- fileLines[!grepl("^__", fileLines)]

    ## Remove quotation marks
    fileLines <- gsub('"', '', fileLines)

    ## Create output filename
    outputFilePath <- sub("\\.txt$", ".clean.txt", inputFilePath)
    writeLines(fileLines, outputFilePath)

    return(outputFilePath)
}

## Filter exons based on quality criteria --------------------------------------
filterExonsByQuality <- function(dexseqObject, minExonLength) {
    ## Filter 1: At least one sample with feature counts >= 50
    keepByMinCount <- rowSums(featureCounts(dexseqObject) >= 50) > 0

    ## Filter 2: All "others" columns must have counts >= 10
    splitColumns <- split(seq_len(ncol(dexseqObject)), colData(dexseqObject)$exon)
    otherCounts <- counts(dexseqObject)[, splitColumns$others]
    keepByOtherCounts <- rowSums(otherCounts >= 10) == ncol(otherCounts)

    ## Filter 3: Exon length >= minimum threshold
    exonLengths <- dexseqObject@rowRanges@ranges@width
    names(exonLengths) <- rownames(mcols(dexseqObject))
    keepByLength <- exonLengths >= minExonLength

    ## Combine all filters
    keepRows <- keepByMinCount & keepByOtherCounts & keepByLength

    return(dexseqObject[keepRows, ])
}

## =============================================================================
## Prepare DEXSeq dataset
## =============================================================================

## Create output directory for top DEU plots -----------------------------------
topDeuPlotDir <- "top_DEU_plots"
dir.create(topDeuPlotDir, showWarnings = FALSE)

## Clean all count files -------------------------------------------------------
countFilePattern <- paste0(samplesheet$sample_id, ".counts.txt")
cleanedCountFiles <- sapply(countFilePattern, cleanDexseqCountFile)

## Create DEXSeq dataset -------------------------------------------------------
dexseqDataset <- DEXSeqDataSetFromHTSeq(
    countfiles = cleanedCountFiles,
    sampleData = samplesheet,
    design = ~ sample + exon + condition:exon,
    flattenedfile = flattenedGffFile
)

## =============================================================================
## Filter exons and run analysis
## =============================================================================

## Apply quality filters -------------------------------------------------------
dexseqDatasetFiltered <- filterExonsByQuality(dexseqDataset, minExonLength)

## Run DEXSeq analysis ---------------------------------------------------------
dexseqResults <- DEXSeq(dexseqDatasetFiltered)

## Remove NA values ------------------------------------------------------------
dexseqResultsClean <- na.omit(dexseqResults)

## =============================================================================
## Identify significant differential exon usage
## =============================================================================

## Apply thresholds from command line ------------------------------------------
significantIndices <- dexseqResultsClean$padj < pvalueThreshold &
                      abs(dexseqResultsClean[, 10]) > logFcThreshold

## Store results in structured list --------------------------------------------
resultsList <- list(
    all = dexseqResultsClean,
    significant = dexseqResultsClean[significantIndices, ]
)

## Prepare data frames for output ----------------------------------------------
resultsAllDataFrame <- data.frame(
    group_id = dexseqResultsClean$groupID,
    feature_id = dexseqResultsClean$featureID,
    log2_fc = dexseqResultsClean[, 10],
    pvalue = dexseqResultsClean$pvalue,
    padj = dexseqResultsClean$padj,
    stringsAsFactors = FALSE
)

resultsSignificantDataFrame <- data.frame(
    group_id = dexseqResultsClean[significantIndices, ]$groupID,
    feature_id = dexseqResultsClean[significantIndices, ]$featureID,
    log2_fc = dexseqResultsClean[significantIndices, 10],
    pvalue = dexseqResultsClean[significantIndices, ]$pvalue,
    padj = dexseqResultsClean[significantIndices, ]$padj,
    stringsAsFactors = FALSE
)

## =============================================================================
## Save results
## =============================================================================

## Save CSV files --------------------------------------------------------------
write.csv(resultsAllDataFrame, "DEXSeq_table.csv", row.names = TRUE)
write.csv(resultsSignificantDataFrame, "significant_DEUs.csv", row.names = TRUE)

## Save significant gene IDs ---------------------------------------------------
saveRDS(
    unique(resultsList$significant$groupID),
    "significant_DEU_gene_ids.rds"
)

## =============================================================================
## Create visualizations
## =============================================================================

## MA plot ---------------------------------------------------------------------
pdf("DEXSeq_MA_plot.pdf", width = 10, height = 8)
plotMA(resultsList$all, cex = 0.8)
dev.off()

## Generate plots for top genes ------------------------------------------------
geneIds <- head(
    unique(resultsList$significant$groupID[order(resultsList$significant$padj)]),
    n = topNGenes
)

for (geneId in geneIds) {
    pdfFilePath <- file.path(topDeuPlotDir, paste0("DEXSeq_", geneId, "_plot.pdf"))
    pdf(pdfFilePath, width = 14, height = 12)

    tryCatch({
        ## Plot 1: Exon usage
        plotDEXSeq(
            resultsList$all,
            geneId,
            legend = TRUE,
            cex.axis = 1.2,
            cex = 1.3,
            lwd = 2,
            FDR = pvalueThreshold,
            displayTranscripts = TRUE,
            main = paste("Gene:", geneId)
        )

        ## Plot 2: Normalized counts
        plotDEXSeq(
            resultsList$all,
            geneId,
            legend = TRUE,
            cex.axis = 1.2,
            cex = 1.3,
            lwd = 2,
            expression = FALSE,
            norCounts = TRUE,
            FDR = pvalueThreshold,
            displayTranscripts = TRUE,
            main = paste("Gene:", geneId)
        )

        ## Plot 3: Splicing
        plotDEXSeq(
            resultsList$all,
            geneId,
            legend = TRUE,
            cex.axis = 1.2,
            cex = 1.3,
            lwd = 2,
            expression = FALSE,
            splicing = TRUE,
            FDR = pvalueThreshold,
            displayTranscripts = TRUE,
            main = paste("Gene:", geneId)
        )
    }, error = function(error) {
        plot(1, 1, type = "n", main = paste("Failed:", geneId))
        text(1, 1, paste("Error:", error$message), col = "red")
    })
    dev.off()
}

## =============================================================================
## Generate HTML report
## =============================================================================

DEXSeqHTML(
    resultsList$all,
    genes = unique(resultsList$significant$groupID),
    FDR = pvalueThreshold,
    color = c("#FF000080", "#0000FF80")
)
