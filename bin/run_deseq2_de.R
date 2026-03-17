#!/usr/bin/env Rscript

## =============================================================================
## Script:  run_deseq2_de.R
## Purpose: Perform differential expression analysis using DESeq2
## Usage:   Rscript run_deseq2_de.R <count_matrix.rds> <samplesheet.csv> <pvalue> <logfc> <top_n>
## Example: Rscript run_deseq2_de.R counts.rds samples.csv 0.05 1.0 100
## =============================================================================

## Load required packages ------------------------------------------------------
library(DESeq2)

## Parse command line arguments ------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
    stop("Usage: run_deseq2_de.R <count_matrix.rds> <samplesheet.csv> <pvalue> <logfc> <top_n>",
         call. = FALSE)
}

countMatrixFile <- args[1]
samplesheetFile <- args[2]
pValueThreshold <- as.numeric(args[3])
logFcThreshold <- as.numeric(args[4])
topNGenes <- as.numeric(args[5])

## =============================================================================
## Load and validate input data
## =============================================================================

## Load count matrix (RDS file)
countData <- readRDS(countMatrixFile)

## Load sample information
samplesheet <- read.csv(samplesheetFile, stringsAsFactors = FALSE)
samplesheet <- samplesheet[, c("sample_id", "condition"), drop = FALSE]
samplesheet$condition <- factor(samplesheet$condition)

## =============================================================================
## Prepare DESeq2 dataset
## =============================================================================

## Handle different input formats (tximport output vs raw count matrix)
if ("counts" %in% names(countData)) {
    ## tximport output format - extract counts matrix
    countData$counts <- countData$counts[, samplesheet$sample_id, drop = FALSE]
    countData$abundance <- countData$abundance[, samplesheet$sample_id, drop = FALSE]
    countMatrix <- round(countData$counts)
} else {
    ## Raw count matrix format
    countMatrix <- countData[, samplesheet$sample_id, drop = FALSE]
}

## Create DESeq2 dataset
deseqDataset <- DESeqDataSetFromMatrix(
    countData = countMatrix,
    colData = samplesheet,
    design = ~ condition
)

## =============================================================================
## Filter lowly expressed genes
## =============================================================================

## Save unfiltered log2 counts
saveRDS(log2(counts(deseqDataset) + 1), file = "unfiltered_counts.rds")

## Filter genes with counts > 10 in at least half of the samples
keepGenes <- rowSums(counts(deseqDataset) > 10) >= (ncol(deseqDataset) / 2)
deseqDataset <- deseqDataset[keepGenes, ]

## Save filtered log2 counts
saveRDS(log2(counts(deseqDataset) + 1), file = "filtered_counts.rds")

## =============================================================================
## Define contrast levels
## =============================================================================

if ("control" %in% levels(samplesheet$condition)) {
    referenceLevel <- "control"
    otherLevel <- setdiff(levels(samplesheet$condition), referenceLevel)[1]
} else {
    referenceLevel <- levels(samplesheet$condition)[1]
    otherLevel <- levels(samplesheet$condition)[2]
}

## Relevel factor
deseqDataset$condition <- relevel(deseqDataset$condition, ref = referenceLevel)

## =============================================================================
## Run DESeq2 analysis with error handling
## =============================================================================

deseqDataset <- tryCatch({
    DESeq(deseqDataset)
}, error = function(e) {
    if (grepl("dispersion estimates", e$message)) {
        deseqDataset <- estimateSizeFactors(deseqDataset)
        deseqDataset <- estimateDispersionsGeneEst(deseqDataset)
        dispersions(deseqDataset) <- mcols(deseqDataset)$dispGeneEst
        deseqDataset <- nbinomWaldTest(deseqDataset)
        return(deseqDataset)
    } else {
        stop(e)
    }
})

## =============================================================================
## Extract normalized counts
## =============================================================================

normalizedCounts <- counts(deseqDataset, normalized = TRUE)
saveRDS(log2(normalizedCounts + 1), file = "normalized_counts.rds")

## =============================================================================
## Extract results
## =============================================================================

resultsObject <- results(deseqDataset, contrast = c("condition", otherLevel, referenceLevel))

## Apply log-fold change shrinkage
resultsShrunk <- tryCatch({
    lfcShrink(deseqDataset, coef = 2, type = "apeglm")
}, error = function(e) {
    return(resultsObject)
})

## =============================================================================
## Variance-stabilizing transformation
## =============================================================================

vstObject <- tryCatch({
    if (nrow(normalizedCounts) < 5000) {
        varianceStabilizingTransformation(deseqDataset, blind = FALSE)
    } else {
        vst(deseqDataset, blind = FALSE)
    }
}, error = function(e) {
    if (grepl("dispersion estimates", e$message)) {
        return(normalizedCounts)
    } else {
        stop(e)
    }
})

## =============================================================================
## Identify differentially expressed genes
## =============================================================================

## Prepare results data frame
resultsDataFrame <- as.data.frame(resultsShrunk)
colnames(resultsDataFrame)[colnames(resultsDataFrame) == "log2FoldChange"] <- "logFC"
colnames(resultsDataFrame)[colnames(resultsDataFrame) == "padj"] <- "adj.P.Val"

## Identify significant genes
significantGenes <- rownames(resultsDataFrame)[
    !is.na(resultsDataFrame$adj.P.Val) &
    resultsDataFrame$adj.P.Val < pValueThreshold &
    abs(resultsDataFrame$logFC) > logFcThreshold
]

## Extract VST values for significant genes
if (inherits(vstObject, "DESeqTransform")) {
    diffExpMatrix <- assay(vstObject)[significantGenes, , drop = FALSE]
} else {
    diffExpMatrix <- vstObject[significantGenes, , drop = FALSE]
}

## Create table of significant genes
diffExpTable <- resultsDataFrame[significantGenes, , drop = FALSE]
diffExpTable <- diffExpTable[order(-abs(diffExpTable$logFC)), ]

## Extract top N genes
topGenes <- head(diffExpTable, n = topNGenes)

## =============================================================================
## Save output files
## =============================================================================

saveRDS(vstObject, file = "vst_transformed.rds")
saveRDS(diffExpMatrix, file = "deseq2_de_genes.rds")
saveRDS(resultsDataFrame, file = "unfiltered_deg_results.rds")
saveRDS(diffExpTable, file = "de_genes_full_table.rds")
saveRDS(topGenes, file = paste0("top_", topNGenes, "_genes.rds"))

write.csv(topGenes, file = paste0("top_", topNGenes, "_genes.csv"), row.names = TRUE)
