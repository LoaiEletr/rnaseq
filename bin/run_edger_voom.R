#!/usr/bin/env Rscript

## =============================================================================
## Script:  run_edger_voom.R
## Purpose: Perform edgeR analysis with voom transformation for RNA-seq data
## Usage:   Rscript run_edger_voom.R <count_matrix.rds> <samplesheet.csv>
## =============================================================================

## Load required packages ------------------------------------------------------
library(edgeR)

## Parse command line arguments ------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Usage: run_edger_voom.R <count_matrix.rds> <samplesheet.csv>",
         call. = FALSE)
}

countMatrixFile <- args[1]
samplesheetFile <- args[2]

## =============================================================================
## Load and validate input data
## =============================================================================

## Load count matrix (RDS file)
countData <- readRDS(countMatrixFile)

## Load sample information
samplesheet <- read.csv(samplesheetFile, stringsAsFactors = FALSE)

## Extract required columns
samplesheet <- samplesheet[, c("sample_id", "condition"), drop = FALSE]
samplesheet$condition <- factor(samplesheet$condition)

## =============================================================================
## Prepare DGEList object
## =============================================================================

## Handle different input formats (tximport output vs raw count matrix)
if ("counts" %in% names(countData)) {
    ## tximport output format - extract counts matrix
    countData$counts <- countData$counts[, samplesheet$sample_id, drop = FALSE]
    countData$abundance <- countData$abundance[, samplesheet$sample_id, drop = FALSE]
    dgeList <- DGEList(countData$counts)
} else {
    ## Raw count matrix format
    countData <- countData[, samplesheet$sample_id, drop = FALSE]
    dgeList <- DGEList(countData)
}

## =============================================================================
## Filter lowly expressed genes
## =============================================================================

## Calculate CPM values
cpmValues <- cpm(dgeList)

## Keep genes with CPM > 10 in at least half of the samples
keepGenes <- rowSums(cpmValues > 10) >= (ncol(cpmValues) / 2)
dgeListFiltered <- dgeList[keepGenes, , keep.lib.sizes = FALSE]

## =============================================================================
## TMM normalization
## =============================================================================

dgeListNormalized <- calcNormFactors(dgeListFiltered, method = "TMM")

## =============================================================================
## Voom transformation with error handling
## =============================================================================

pdf("voom_plot.pdf")
if (nrow(dgeListFiltered$counts) >= 2) {
    ## Create design matrix
    group <- factor(samplesheet$condition)
    designMatrix <- model.matrix(~0 + group)
    colnames(designMatrix) <- levels(group)

    ## Apply voom transformation
    voomObject <- voom(dgeListNormalized, designMatrix, plot = TRUE)
} else {
    ## Error plot if insufficient genes
    plot(1, 1,
         type = "n",
         main = "Error: Insufficient genes for voom transformation",
         xlab = "", ylab = "",
         xaxt = "n", yaxt = "n")
    text(1, 1,
         paste("Genes after filtering:", nrow(dgeListFiltered$counts)),
         cex = 1.5)
    voomObject <- NULL
}
dev.off()

## =============================================================================
## Save output files
## =============================================================================

## Save log2 CPM matrices at different stages
saveRDS(cpm(dgeList, log = TRUE), file = "unfiltered_log2cpm.rds")
saveRDS(cpm(dgeListFiltered, log = TRUE), file = "filtered_log2cpm.rds")
saveRDS(cpm(dgeListNormalized, log = TRUE), file = "normalized_log2cpm.rds")

## Save voom object if successful
if (!is.null(voomObject)) {
    saveRDS(voomObject, file = "voom_transformed.rds")
    cat("Voom transformation successful\n")
} else {
    cat("WARNING: Voom transformation failed - insufficient genes\n")
}
