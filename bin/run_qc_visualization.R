#!/usr/bin/env Rscript

## =============================================================================
## Script:  run_qc_visualization.R
## Purpose: Generate all QC plots (MDS, PCA, Violin) in one module
## Usage:   Rscript run_qc_visualization.R <voom_dgelist> <samplesheet_csv> <log2_unfiltered> <log2_filtered> <log2_normalized> <vst_object>
## Example: Rscript run_qc_visualization.R voom.rds samples.csv unfiltered.rds filtered.rds normalized.rds vst.rds
## =============================================================================

library(edgeR)
library(ggplot2)
library(DESeq2)
library(tibble)
library(tidyr)
library(dplyr)
library(cowplot)
library(MASS)

## =============================================================================
## Helper Functions
## =============================================================================

## Create MDS plot
createMdsPlot <- function(voomDgeList, samplesheet) {
    if (is.null(voomDgeList)) {
        plot.new()
        text(
            x = 0.5,
            y = 0.5,
            labels = "MDS plot not generated because the voom object is NULL.",
            cex = 1.2
        )
        return(NULL)
    }

    plotMDS(
        voomDgeList,
        labels = samplesheet$condition,
        col    = as.numeric(samplesheet$condition),
        main   = "MDS Plot - Samples by Group"
    )
}

## Create PCA plot with VST fallback
createPcaPlot <- function(vstObject, samplesheet, normalizedCounts = NULL) {

    samplesheet$condition <- factor(samplesheet$condition)

    ## Try to use VST object first
    if (inherits(vstObject, "DESeqTransform")) {
        cat("Using VST object for PCA plot\n")

        ## Create PCA plot using DESeq2's built-in function
        plotPCA(vstObject, intgroup = "condition")

    } else if (!is.null(normalizedCounts)) {
        cat("Using normalized counts for PCA plot\n")

        ## Convert to matrix for normalized data
        countMatrix <- normalizedCounts

        ## Perform standard PCA
        tryCatch({
            pcaRes <- prcomp(t(countMatrix), scale. = FALSE, center = TRUE)
            pcaDf <- as.data.frame(pcaRes$x)

            if (ncol(pcaDf) < 2) {
                plot.new()
                text(0.5, 0.5, "Less than 2 principal components", cex = 1.2)
            } else {
                varPer <- round(100 * pcaRes$sdev^2 / sum(pcaRes$sdev^2), 1)

                pcaDf$sample_id <- samplesheet$sample_id
                pcaDf$group <- samplesheet$condition

                p <- ggplot(pcaDf, aes(x = PC1, y = PC2, color = group)) +
                    geom_point(size = 4) +
                    stat_ellipse() +
                    xlab(paste0("PC1 (", varPer[1], "%)")) +
                    ylab(paste0("PC2 (", varPer[2], "%)")) +
                    ggtitle("PCA Plot") +
                    theme_bw() +
                    coord_fixed()

                if (nrow(pcaDf) <= 30) {
                    p <- p + geom_text(aes(label = sample_id),
                                       hjust = 0.5, vjust = -0.5, size = 3,
                                       show.legend = FALSE)
                }

                print(p)
            }
        }, error = function(e) {
            plot.new()
            text(0.5, 0.5, paste("PCA failed:", e$message), cex = 1.2)
        })
    } else {
        plot.new()
        text(0.5, 0.5, "No valid data for PCA plot", cex = 1.2)
    }
}

## Process log2 counts data
readAndProcessLog2Counts <- function(filePath, description) {
    ## Read the data
    log2Data <- readRDS(filePath)

    ## Convert to tibble and pivot
    log2Df <- as_tibble(log2Data, rownames = "geneId")
    log2DfPivot <- pivot_longer(log2Df,
                                cols = -geneId,
                                names_to = "samples",
                                values_to = "expression")

    return(list(data = log2DfPivot, description = description))
}

## Create violin plot
createViolinPlot <- function(data, description) {
    ggplot(data) +
        aes(x = samples, y = expression, fill = samples) +
        geom_violin(trim = FALSE, show.legend = FALSE) +
        stat_summary(fun = "median",
                     geom = "point",
                     shape = 95,
                     size = 10,
                     color = "black",
                     show.legend = FALSE) +
        labs(y = "log2 expression",
             x = "sample",
             title = "Log2 Counts",
             subtitle = description,
             caption = paste0("produced on ", Sys.time())) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

## =============================================================================
## Main Execution
## =============================================================================

## Parse command line arguments ------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
    cat("Usage: run_qc_visualization.R <voom_dgelist> <samplesheet_csv> <log2_unfiltered> <log2_filtered> <log2_normalized> <vst_object>\n")
    cat("Arguments received:", length(args), "\n")
    stop("Insufficient arguments", call. = FALSE)
}

voomDgeListFile    <- args[1]
samplesheetFile    <- args[2]
log2UnfilteredFile <- args[3]
log2FilteredFile   <- args[4]
log2NormalizedFile <- args[5]
vstObjectFile      <- args[6]

## Load data -------------------------------------------------------------------

voomDgeList <- if (voomDgeListFile != "null") readRDS(voomDgeListFile) else NULL
samplesheet <- read.csv(samplesheetFile)
samplesheet <- samplesheet[, c("sample_id", "condition"), drop = FALSE]

## Read all log2 datasets
unfilteredData <- readAndProcessLog2Counts(log2UnfilteredFile,
                                           "unfiltered, non-normalized")
filteredData <- readAndProcessLog2Counts(log2FilteredFile,
                                         "filtered, non-normalized")
normalizedData <- readAndProcessLog2Counts(log2NormalizedFile,
                                           "filtered, normalized")

## Load VST object
vstObject <- readRDS(vstObjectFile)

## =============================================================================
## Generate all QC plots
## =============================================================================

## Create output directory
qcDir <- "qc_plots"
dir.create(qcDir, showWarnings = FALSE)

## 1. MDS Plot
pdf(file.path(qcDir, "mds_plot.pdf"))
createMdsPlot(voomDgeList, samplesheet)
dev.off()

## 2. PCA Plot (using VST object)
pdf(file.path(qcDir, "pca_plot.pdf"))
createPcaPlot(vstObject, samplesheet, readRDS(log2NormalizedFile))
dev.off()

## 3. Violin Plots
p1 <- createViolinPlot(unfilteredData$data, unfilteredData$description)
p2 <- createViolinPlot(filteredData$data, filteredData$description)
p3 <- createViolinPlot(normalizedData$data, normalizedData$description)

## Combine plots
finalPlot <- plot_grid(p1, p2, p3,
                      labels = c('A', 'B', 'C'),
                      label_size = 12,
                      ncol = 1)

## Save the plot
ggsave(file.path(qcDir, "violin_plots.pdf"), finalPlot, width = 12, height = 15)
