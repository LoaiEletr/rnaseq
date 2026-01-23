#!/usr/bin/env Rscript

## =============================================================================
## Script:  run_limma_de.R
## Purpose: Perform differential expression analysis using limma
## Usage:   Rscript run_limma_de.R <voom_object.rds> <samplesheet.csv> <pvalue> <logfc> <top_n>
## Example: Rscript run_limma_de.R voom_transformed.rds samples.csv 0.05 1.0 100
## =============================================================================

## Load required packages ------------------------------------------------------
library(limma)

## Parse command line arguments ------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
    stop("Usage: run_limma_de.R <voom_object.rds> <samplesheet.csv> <pvalue> <logfc> <top_n>",
         call. = FALSE)
}

voomObjectFile <- args[1]
samplesheetFile <- args[2]
pValueThreshold <- as.numeric(args[3])
logFcThreshold <- as.numeric(args[4])
topNGenes <- as.numeric(args[5])

## =============================================================================
## Load and validate input data
## =============================================================================

## Load voom-transformed object
voomObject <- readRDS(voomObjectFile)

## Load sample information
samplesheet <- read.csv(samplesheetFile, stringsAsFactors = FALSE)
samplesheet <- samplesheet[, c("sample_id", "condition"), drop = FALSE]

## =============================================================================
## Case 1: Voom object is NULL - skip analysis
## =============================================================================
if (is.null(voomObject)) {
    reasonMessage <- paste(
        "Top genes not generated because the voom object is NULL.",
        "Differential expression analysis was skipped."
    )

    ## Save NULL outputs (expected by workflow)
    saveRDS(NULL, file = "ebfit.rds")
    saveRDS(NULL, file = "unfiltered_deg_results.rds")
    saveRDS(NULL, file = "de_expression_values.rds")
    saveRDS(NULL, file = "de_genes_full_table.rds")
    saveRDS(NULL, file = paste0("top_", topNGenes, "_genes.rds"))

    ## Save reason as CSV (not empty)
    reasonDataFrame <- data.frame(
        reason = reasonMessage,
        stringsAsFactors = FALSE
    )

    write.csv(
        reasonDataFrame,
        file = paste0("top_", topNGenes, "_genes.csv"),
        row.names = FALSE
    )

    quit(save = "no", status = 0)
}

## =============================================================================
## Case 2: Normal differential expression analysis
## =============================================================================

## Prepare design matrix -------------------------------------------------------
group <- factor(samplesheet$condition)
designMatrix <- model.matrix(~0 + group)
colnames(designMatrix) <- levels(group)

## Fit linear model ------------------------------------------------------------
fit <- lmFit(voomObject, designMatrix)

## Define contrast -------------------------------------------------------------
if ("control" %in% levels(group)) {
    referenceLevel <- "control"
    otherLevel <- setdiff(levels(group), referenceLevel)[1]
} else {
    referenceLevel <- levels(group)[1]
    otherLevel <- levels(group)[2]
}

contrastMatrix <- makeContrasts(
    contrast = paste0(otherLevel, "-", referenceLevel),
    levels = designMatrix
)

## Apply contrast and empirical Bayes ------------------------------------------
contrastFit <- contrasts.fit(fit, contrastMatrix)
ebayesFit <- eBayes(contrastFit)

## Identify differentially expressed genes -------------------------------------
decideTestsResult <- decideTests(
    ebayesFit,
    method = "global",
    adjust.method = "BH",
    p.value = pValueThreshold,
    lfc = logFcThreshold
)

## Extract differentially expressed genes --------------------------------------
diffExpGenes <- voomObject$E[decideTestsResult[, 1] != 0, , drop = FALSE]

## Get full results table ------------------------------------------------------
fullResults <- topTable(
    ebayesFit,
    adjust = "BH",
    coef = 1,
    number = Inf,
    sort.by = "none"
)

## Filter to significant genes -------------------------------------------------
diffExpTable <- fullResults[
    rownames(fullResults) %in% rownames(diffExpGenes), , drop = FALSE
]

## Get top N genes by logFC ----------------------------------------------------
topGenes <- topTable(
    ebayesFit,
    adjust = "BH",
    coef = 1,
    n = topNGenes,
    sort.by = "logFC"
)

## =============================================================================
## Save output files
## =============================================================================

saveRDS(ebayesFit, file = "ebfit.rds")
saveRDS(fullResults, file = "unfiltered_deg_results.rds")
saveRDS(diffExpGenes, file = "de_expression_values.rds")
saveRDS(diffExpTable, file = "de_genes_full_table.rds")
saveRDS(topGenes, file = paste0("top_", topNGenes, "_genes.rds"))

write.csv(
    topGenes,
    file = paste0("top_", topNGenes, "_genes.csv"),
    row.names = TRUE
)
