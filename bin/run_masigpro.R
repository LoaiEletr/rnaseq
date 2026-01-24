#!/usr/bin/env Rscript

## =============================================================================
## Script:  run_masigpro.R
## Purpose: Perform time-course differential expression analysis using maSigPro
## Usage:   Rscript run_masigpro.R <count_matrix.rds> <samplesheet.csv> <pvalue> <rsq> <cluster_method>
## Example: Rscript run_masigpro.R counts.rds samples.csv 0.05 0.7 hclust
## Cluster methods: hclust, Mclust, kmeans
## =============================================================================

## Load required packages ------------------------------------------------------
library(maSigPro)
library(edgeR)
library(dplyr)
library(tibble)
library(ggplot2)
library(reshape2)
library(patchwork)
library(mclust)

## Set random seed for reproducibility ----------------------------------------
set.seed(1234)

## Parse command line arguments ------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
    stop("Usage: run_masigpro.R <count_matrix.rds> <samplesheet.csv> <pvalue> <rsq> <cluster_method>",
         call. = FALSE)
}

countMatrixFile <- args[1]
samplesheetFile <- args[2]
pValueThreshold <- as.numeric(args[3])
rsqThreshold <- as.numeric(args[4])
clusterMethod <- args[5]

## =============================================================================
## Validate cluster method
## =============================================================================
validClusterMethods <- c("hclust", "Mclust", "kmeans")
if (!clusterMethod %in% validClusterMethods) {
    stop("Invalid cluster_method: ", clusterMethod,
         ". Must be one of: hclust, Mclust, kmeans.", call. = FALSE)
}

## =============================================================================
## Load and validate input data
## =============================================================================

## Load count matrix (RDS file)
countData <- readRDS(countMatrixFile)

## Load sample information
samplesheet <- read.csv(samplesheetFile, stringsAsFactors = FALSE)

## Validate required columns
requiredCols <- c("sample_id", "condition", "timepoint")
if (!all(requiredCols %in% colnames(samplesheet))) {
    stop("Samplesheet must contain columns: ", paste(requiredCols, collapse = ", "),
         call. = FALSE)
}

## Extract required columns
samplesheet <- samplesheet[, requiredCols, drop = FALSE]
samplesheet$timepoint <- as.numeric(samplesheet$timepoint)

## =============================================================================
## Prepare count matrix
## =============================================================================

## Handle different input formats (tximport output vs raw count matrix)
if ("counts" %in% names(countData)) {
    ## tximport output format - extract counts matrix
    countData$counts <- countData$counts[, samplesheet$sample_id, drop = FALSE]
    countData$abundance <- countData$abundance[, samplesheet$sample_id, drop = FALSE]
    countMatrix <- countData$counts
} else {
    ## Raw count matrix format
    countMatrix <- countData[, samplesheet$sample_id, drop = FALSE]
}

## =============================================================================
## Create experimental design matrix for maSigPro
## =============================================================================

designMatrix <- samplesheet %>%
    group_by(condition) %>%
    mutate(Replicate = ave(seq_len(n()), condition, FUN = seq_along)) %>%
    ungroup() %>%
    mutate(Disease = ifelse(condition != "control", 1, 0),
           Control = ifelse(condition == "control", 1, 0)) %>%
    select(-condition) %>%
    column_to_rownames("sample_id") %>%
    select(timepoint, Replicate, Control, Disease) %>%
    as.data.frame()

## Determine polynomial degree (timepoints - 1)
polyDegree <- length(unique(designMatrix$timepoint)) - 1

## =============================================================================
## Prepare maSigPro design
## =============================================================================

maSigProDesign <- make.design.matrix(
    designMatrix,
    degree = polyDegree,
    time.col = 1,
    repl.col = 2,
    group.cols = 3:ncol(designMatrix)
)

## =============================================================================
## EdgeR preprocessing: Filtering and normalization
## =============================================================================

## Create DGEList object
dgeList <- DGEList(countMatrix)

## Filter lowly expressed genes (CPM > 10 in at least half of samples)
keepGenes <- rowSums(cpm(dgeList) >= 10) >= (ncol(countMatrix) / 2)
dgeListFiltered <- dgeList[keepGenes, ]

## TMM normalization
dgeListNormalized <- calcNormFactors(dgeListFiltered, method = "TMM")

## =============================================================================
## maSigPro differential expression analysis
## =============================================================================

## Step 1: Identify significant genes using p.vector
sigProFit <- p.vector(
    cpm(dgeListNormalized),
    design = maSigProDesign,
    Q = pValueThreshold,
    MT.adjust = "BH",
    min.obs = 20
)

## Step 2: Model selection using T.fit
sigProTStep <- T.fit(
    data = sigProFit,
    step.method = "backward",
    alfa = pValueThreshold
)

## Step 3: Extract significant genes based on R-squared threshold
significantGenes <- get.siggenes(
    sigProTStep,
    rsq = rsqThreshold,
    vars = "groups"
)

sigGeneNames <- rownames(significantGenes$sig.genes$DiseasevsControl$sig.profiles)

## =============================================================================
## Clustering of significant genes
## =============================================================================

## Create output directory
outputDir <- "masigpro_results"
dir.create(outputDir, showWarnings = FALSE)

## Prepare clustering parameters
useMclust <- ifelse(clusterMethod %in% c("hclust", "kmeans"), FALSE, TRUE)

clusterParams <- list(
    data = significantGenes$sig.genes$DiseasevsControl,
    edesign = maSigProDesign$edesign,
    cluster.method = clusterMethod,
    cluster.data = 1,
    show.lines = FALSE,
    show.fit = TRUE,
    summary.mode = "median",
    min.obs = 7,
    dis = maSigProDesign$dis,
    k.mclust = useMclust,
    groups.vector = maSigProDesign$groups.vector,
    alfa = pValueThreshold,
    color.mode = "rainbow",
    lwd = 2,
    pch = 19,
    cex = 1.3,
    legend = TRUE,
    x.legend = "topright",
    cex.legend = 0.8
)

## Add k parameter for non-Mclust methods
if (!useMclust) {
    clusterParams$k <- 9
}

## Perform clustering
par(mfrow = c(3, 3), mar = c(5, 5, 4, 3), oma = c(2, 2, 2, 2))
geneClusters <- do.call(see.genes, clusterParams)

## =============================================================================
## Prepare clustering results for visualization
## =============================================================================

## Extract cluster assignments
clusterAssignments <- data.frame(
    gene_id = names(geneClusters$cut),
    cluster = as.integer(geneClusters$cut),
    stringsAsFactors = FALSE
)

## Get expression data for significant genes
sigGeneExpression <- cpm(dgeListNormalized)[sigGeneNames, , drop = FALSE]

## Reshape expression data for plotting
expressionMelted <- reshape2::melt(
    as.matrix(sigGeneExpression),
    varnames = c("gene_id", "sample_id"),
    value.name = "cpm"
)

## Add metadata to melted data
expressionMelted$condition <- samplesheet$condition[match(expressionMelted$sample_id, samplesheet$sample_id)]
expressionMelted$timepoint <- samplesheet$timepoint[match(expressionMelted$sample_id, samplesheet$sample_id)]
expressionMelted$cluster <- geneClusters$cut[expressionMelted$gene_id]

## Calculate summary statistics per cluster
clusterSummary <- expressionMelted %>%
    group_by(cluster, condition, timepoint) %>%
    summarise(median_cpm = median(cpm), .groups = "drop") %>%
    arrange(cluster, condition, timepoint)

clusterGeneCounts <- expressionMelted %>%
    group_by(cluster) %>%
    summarise(n_genes = n_distinct(gene_id), .groups = "drop")

## =============================================================================
## Create cluster visualization plots
## =============================================================================

## Define color scheme
conditions <- unique(samplesheet$condition)
nonControlConditions <- setdiff(conditions, "control")
conditionColors <- c(control = "#1F78B4", setNames("#E31A1C", nonControlConditions))

## Function to plot individual cluster
plotCluster <- function(summaryData, countData, clusterNum, showLegend = TRUE) {
    nGenes <- countData$n_genes[countData$cluster == clusterNum]
    plotData <- summaryData[summaryData$cluster == clusterNum, ]

    ggplot(plotData, aes(x = timepoint, y = median_cpm, color = condition, group = condition)) +
        geom_line(linewidth = 1) +
        geom_point(size = 2.5) +
        scale_color_manual(values = conditionColors, name = "Condition") +
        scale_x_continuous(breaks = unique(plotData$timepoint)) +
        labs(title = paste0("Cluster ", clusterNum, " (n=", nGenes, ")"),
             x = "Timepoint", y = "Median CPM") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = ifelse(showLegend, "right", "none"))
}

## Create plots for all clusters
allClusters <- sort(unique(clusterSummary$cluster))
clusterPlots <- lapply(allClusters, function(cl) {
    plotCluster(clusterSummary, clusterGeneCounts, cl, showLegend = TRUE)
})

## Combine all plots
combinedPlot <- wrap_plots(clusterPlots, ncol = 3) +
    plot_layout(guides = "collect")

## Save combined plot
ggsave(file.path(outputDir, "masigpro_clusterplots.pdf"),
       combinedPlot, width = 15, height = 12, dpi = 300)

## =============================================================================
## Save output files
## =============================================================================

## Summary statistics
summaryStats <- data.frame(
    total_genes_tested = nrow(dgeListNormalized),
    significant_genes = length(sigGeneNames),
    num_clusters = length(clusterGeneCounts$cluster),
    polynomial_degree = polyDegree,
    q_threshold = pValueThreshold,
    rsq_threshold = rsqThreshold,
    cluster_method = clusterMethod
)

## Save CSV files
write.csv(clusterAssignments, file.path(outputDir, "cluster_assignments.csv"), row.names = FALSE)
write.csv(clusterSummary, file.path(outputDir, "cluster_median_expression.csv"), row.names = FALSE)
write.csv(summaryStats, file.path(outputDir, "summary_statistics.csv"), row.names = FALSE)

## Save RDS files
saveRDS(sigGeneExpression, file.path(outputDir, "expression_matrix_siggenes.rds"))
saveRDS(geneClusters, file.path(outputDir, "clusters_list.rds"))
saveRDS(cpm(dgeList, log = TRUE), file.path(outputDir, "unfiltered_log2cpm.rds"))
saveRDS(cpm(dgeListFiltered, log = TRUE), file.path(outputDir, "filtered_log2cpm.rds"))
saveRDS(cpm(dgeListNormalized, log = TRUE), file.path(outputDir, "normalized_log2cpm.rds"))
