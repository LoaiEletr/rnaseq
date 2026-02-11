#!/usr/bin/env Rscript

## =============================================================================
## Script:  run_wgcna_analysis.R
## Purpose: Weighted Gene Co-expression Network Analysis (WGCNA)
## Usage:   Rscript run_wgcna_analysis.R <counts_file> <samplesheet_file> <min_gene_significance>
##          <min_module_membership> <top_n_genes> <correlation_threshold> <pvalue_threshold>
##          <max_block_size> <tom_type> <network_type> <deep_split> <reassign_threshold>
##          <merge_cut_height> <min_module_size> <scale_free_r2_threshold>
## Example: Rscript run_wgcna_analysis.R counts.rds samples.csv 0.2 0.5 50 0.5 0.05 5000 signed
##          signed 2 0.05 0.25 30 0.8
## =============================================================================

## Load required packages ------------------------------------------------------
library(DESeq2)
library(WGCNA)
library(pheatmap)

## Enable WGCNA multi-threading ------------------------------------------------
enableWGCNAThreads()

## =============================================================================
## Parse command line arguments
## =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 15) {
    stop(
        "Usage: run_wgcna_analysis.R <counts_file> <samplesheet_file> <min_gene_significance> ",
        "<min_module_membership> <top_n_genes> <correlation_threshold> <pvalue_threshold> ",
        "<max_block_size> <tom_type> <network_type> <deep_split> <reassign_threshold> ",
        "<merge_cut_height> <min_module_size> <scale_free_r2_threshold>",
        call. = FALSE
    )
}

## Extract parameters with descriptive names -----------------------------------
countsFile <- args[1]
samplesheetFile <- args[2]
minGeneSignificance <- as.numeric(args[3])
minModuleMembership <- as.numeric(args[4])
topNGenes <- as.numeric(args[5])
correlationThreshold <- as.numeric(args[6])
pvalueThreshold <- as.numeric(args[7])
maxBlockSize <- as.numeric(args[8])
tomType <- args[9]
networkType <- args[10]
deepSplit <- as.numeric(args[11])
reassignThreshold <- as.numeric(args[12])
mergeCutHeight <- as.numeric(args[13])
minModuleSize <- as.numeric(args[14])
rSquaredThreshold = as.numeric(args[15])

## =============================================================================
## Load and prepare data
## =============================================================================

## Load sample information -----------------------------------------------------
samplesheet <- read.csv(samplesheetFile, stringsAsFactors = FALSE)
samplesheet$condition <- factor(samplesheet$condition)

## Subset samplesheet to required columns --------------------------------------
if (is.null(samplesheet$timepoint)) {
    samplesheet <- samplesheet[, c("sample_id", "condition"), drop = FALSE]
} else {
    samplesheet <- samplesheet[, c("sample_id", "condition", "timepoint"), drop = FALSE]
}

## Load count matrix -----------------------------------------------------------
countData <- readRDS(countsFile)

## Extract counts from tximport object if present ------------------------------
if ("counts" %in% names(countData)) {
    countData$counts <- countData$counts[, samplesheet$sample_id, drop = FALSE]
    countData$abundance <- countData$abundance[, samplesheet$sample_id, drop = FALSE]
    countMatrix <- round(countData$counts)
} else {
    countData <- countData[, samplesheet$sample_id, drop = FALSE]
    countMatrix <- countData
}

## Quality control: remove outlier genes ---------------------------------------
geneSampleGoodness <- goodSamplesGenes(t(countMatrix))

if (!geneSampleGoodness$allOK) {
    nGenesRemoved <- sum(!geneSampleGoodness$goodGenes)
    countMatrix <- countMatrix[geneSampleGoodness$goodGenes == TRUE, ]
    cat("Removed", nGenesRemoved, "outlier genes\n")
}

## =============================================================================
## DESeq2 normalization
## =============================================================================

cat("\n2. Normalizing with DESeq2...\n")

## Create DESeq2 dataset -------------------------------------------------------
if (is.null(samplesheet$timepoint)) {
    deseqDataset <- DESeqDataSetFromMatrix(
        countMatrix,
        colData = samplesheet,
        design = ~ condition
    )
} else {
    samplesheet$timepoint <- factor(samplesheet$timepoint)
    deseqDataset <- DESeqDataSetFromMatrix(
        countMatrix,
        colData = samplesheet,
        design = ~ timepoint + condition + timepoint:condition
    )
}

## Set reference level for condition -------------------------------------------
if ("control" %in% levels(samplesheet$condition)) {
    referenceLevel <- "control"
} else {
    referenceLevel <- levels(samplesheet$condition)[1]
}

deseqDataset$condition <- relevel(deseqDataset$condition, ref = referenceLevel)

## Filter lowly expressed genes ------------------------------------------------
minCountThreshold <- 15
minSampleFraction <- 0.75
minSamplesRequired <- round(ncol(deseqDataset) * minSampleFraction)

keepGenes <- rowSums(counts(deseqDataset) >= minCountThreshold) >= minSamplesRequired
deseqDatasetFiltered <- deseqDataset[keepGenes, ]

## Create output directory -----------------------------------------------------
outputDir <- "WGCNA"
dir.create(outputDir, showWarnings = FALSE)

## Save filtered and unfiltered counts ----------------------------------------
saveRDS(log2(counts(deseqDataset) + 1), file.path(outputDir, "unfiltered_counts.rds"))
saveRDS(log2(counts(deseqDatasetFiltered) + 1), file.path(outputDir, "filtered_counts.rds"))

## Run DESeq2 with error handling ----------------------------------------------
deseqResult <- tryCatch(
    {
        DESeq(deseqDatasetFiltered)
    },
    error = function(error) {
        if (grepl("dispersion estimates", error$message)) {
            message("Dispersion fit failed — using gene-wise dispersion estimates.")
            deseqDatasetFiltered <- estimateSizeFactors(deseqDatasetFiltered)
            deseqDatasetFiltered <- estimateDispersionsGeneEst(deseqDatasetFiltered)
            dispersions(deseqDatasetFiltered) <- mcols(deseqDatasetFiltered)$dispGeneEst
            deseqDatasetFiltered <- nbinomWaldTest(deseqDatasetFiltered)
            return(deseqDatasetFiltered)
        } else {
            stop(error)
        }
    }
)

## Variance stabilization transformation ---------------------------------------
vstResult <- tryCatch(
    {
        if (nrow(deseqResult) < 5000) {
            varianceStabilizingTransformation(deseqResult)
        } else {
            vst(deseqResult)
        }
    },
    error = function(error) {
        errorMessage <- "VST failed — using normalized counts as fallback."
        writeLines(errorMessage, file.path(outputDir, "error_log.txt"))
        stop(error)
    }
)

## Save normalized expression matrix -------------------------------------------
normalizedCounts <- assay(vstResult)
normalizedCounts <- t(normalizedCounts)

saveRDS(vstResult, file.path(outputDir, "vsd_object.rds"))
saveRDS(
    log2(counts(deseqResult, normalized = TRUE) + 1),
    file.path(outputDir, "normalized_counts.rds")
)

## =============================================================================
## Network construction
## =============================================================================

## Soft threshold selection ----------------------------------------------------
powerValues <- c(1:10, seq(from = 12, to = 20, by = 2))
softThresholdResult <- pickSoftThreshold(
    normalizedCounts,
    powerVector = powerValues,
    networkType = "signed",
    verbose = 5
)

## Function to select optimal power --------------------------------------------
selectOptimalPower <- function(fitIndices, rSquaredThreshold = 0.80) {

    for (i in seq_len(nrow(fitIndices))) {
        if (fitIndices$SFT.R.sq[i] >= rSquaredThreshold) {
            selectedPower <- fitIndices$Power[i]
            return(selectedPower)
        }
    }
    selectedPower <- fitIndices$Power[which.max(fitIndices$SFT.R.sq)]
    return(selectedPower)
}

softPower <- selectOptimalPower(softThresholdResult$fitIndices, rSquaredThreshold)

## Plot soft threshold selection -----------------------------------------------
pdf(file.path(outputDir, "soft_threshold_selection.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))

## Scale independence plot
plot(
    softThresholdResult$fitIndices[, 1],
    softThresholdResult$fitIndices[, 2],
    xlab = "Soft Threshold (power)",
    ylab = "Scale Free Topology Model Fit",
    main = "Scale independence",
    type = "n"
)
text(
    softThresholdResult$fitIndices[, 1],
    softThresholdResult$fitIndices[, 2],
    labels = powerValues,
    col = "red"
)
abline(h = 0.8, col = "red")

## Mean connectivity plot
plot(
    softThresholdResult$fitIndices[, 1],
    softThresholdResult$fitIndices[, 5],
    xlab = "Soft Threshold (power)",
    ylab = "Mean Connectivity",
    main = "Mean connectivity",
    type = "n"
)
text(
    softThresholdResult$fitIndices[, 1],
    softThresholdResult$fitIndices[, 5],
    labels = powerValues,
    col = "red"
)
dev.off()

## Network construction with user parameters -----------------------------------
blockwiseNetwork <- blockwiseModules(
    normalizedCounts,
    maxBlockSize = maxBlockSize,
    TOMType = tomType,
    networkType = networkType,
    deepSplit = deepSplit,
    reassignThreshold = reassignThreshold,
    power = softPower,
    mergeCutHeight = mergeCutHeight,
    minModuleSize = minModuleSize,
    numericLabels = FALSE,
    randomSeed = 1234,
    verbose = 3
)

## Plot dendrogram -------------------------------------------------------------
pdf(file.path(outputDir, "module_dendrogram.pdf"), width = 15, height = 8)
plotDendroAndColors(
    blockwiseNetwork$dendrograms[[1]],
    cbind(blockwiseNetwork$unmergedColors, blockwiseNetwork$colors),
    c("Unmerged", "Merged"),
    dendroLabels = FALSE,
    addGuide = TRUE,
    hang = 0.03,
    main = "Gene dendrogram and module colors"
)
dev.off()

## =============================================================================
## Module-trait association
## =============================================================================

cat("\n4. Analyzing module-trait associations...\n")

## Function to create trait matrix from samplesheet ----------------------------
createTimeTraitsCompact <- function(samplesheetData, conditionColumn, timepointColumn) {

    ## Validate inputs
    if (!conditionColumn %in% colnames(samplesheetData)) {
        stop("Column ", conditionColumn, " not found in samplesheet")
    }

    ## Check if timepoint column exists
    if (!timepointColumn %in% colnames(samplesheetData)) {
        ## If no timepoint column, create simple condition columns
        conditions <- unique(samplesheetData[[conditionColumn]])
        traitMatrix <- data.frame(row.names = rownames(samplesheetData))

        for (condition in conditions) {
            if (condition != "control") {
                traitMatrix[[condition]] <- ifelse(
                    samplesheetData[[conditionColumn]] == condition, 1, 0
                )
            }
        }
        return(traitMatrix)
    }

    ## Initialize empty data frame
    traitMatrix <- data.frame(row.names = rownames(samplesheetData))

    ## Get unique combinations of condition and timepoint
    samplesheetData$combination <- paste(
        samplesheetData[[conditionColumn]],
        samplesheetData[[timepointColumn]],
        sep = "_"
    )
    combinations <- unique(samplesheetData$combination)

    ## Create column for each unique combination
    for (combination in combinations) {
        traitMatrix[[combination]] <- ifelse(samplesheetData$combination == combination, 1, 0)
    }

    return(traitMatrix)
}

## Create trait matrix ---------------------------------------------------------
timepointColumn <- if (is.null(samplesheet$timepoint)) "NULL" else "timepoint"

traitMatrix <- createTimeTraitsCompact(
    samplesheetData = samplesheet,
    conditionColumn = "condition",
    timepointColumn = timepointColumn
)

## Get module eigengenes -------------------------------------------------------
moduleEigengenes <- blockwiseNetwork$MEs
moduleEigengenes <- orderMEs(moduleEigengenes)

## Calculate correlations ------------------------------------------------------
nSamples <- nrow(normalizedCounts)
moduleTraitCorrelation <- cor(moduleEigengenes, traitMatrix, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCorrelation, nSamples)

## Module-trait correlation heatmap --------------------------------------------
pdf(file.path(outputDir, "module_trait_correlation_heatmap.pdf"), width = 10, height = 12)

## Create text matrix with significance stars
textMatrix <- paste0(
    signif(moduleTraitCorrelation, 2), "\n",
    ifelse(moduleTraitPvalue < 0.001, "***",
           ifelse(moduleTraitPvalue < 0.01, "**",
                  ifelse(moduleTraitPvalue < 0.05, "*", "")))
)

par(mar = c(8, 10, 3, 3))
labeledHeatmap(
    Matrix = moduleTraitCorrelation,
    xLabels = colnames(moduleTraitCorrelation),
    yLabels = rownames(moduleTraitCorrelation),
    ySymbols = gsub("ME", "", rownames(moduleTraitCorrelation)),
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 0.9,
    zlim = c(-1, 1),
    main = "Module-Trait Relationships",
    cex.main = 1.2
)
dev.off()

rownames(traitMatrix) <- colnames(t(moduleEigengenes))

## Module eigengenes heatmap ---------------------------------------------------
pdf(file.path(outputDir, "module_eigengenes_heatmap.pdf"), width = 12, height = 10)
pheatmap(
    t(moduleEigengenes),
    scale = "row",
    clustering_distance_rows = "euclidean",
    clustering_method = "average",
    annotation_col = traitMatrix,
    show_colnames = FALSE,
    main = "Module Eigengenes Across Samples",
    color = colorRampPalette(c("blue", "white", "red"))(50)
)
dev.off()

## =============================================================================
## Process significant modules for each trait
## =============================================================================

cat("\n7. Processing significant modules for each trait...\n")

## Function to process significant modules -------------------------------------
processSignificantModules <- function(traitColumnIndex, traitName = NULL, nTopGenes = topNGenes) {
    ## If trait_name not provided, use column name
    if (is.null(traitName)) {
        traitName <- colnames(traitMatrix)[traitColumnIndex]
    }

    cat("\n=== Processing trait:", traitName, "===\n")

    ## Identify significant modules for this trait -----------------------------
    significantModules <- rownames(moduleTraitCorrelation)[
        abs(moduleTraitCorrelation[, traitColumnIndex]) > correlationThreshold &
        moduleTraitPvalue[, traitColumnIndex] < pvalueThreshold
    ]

    cat("Found", length(significantModules), "significant modules for", traitName, "\n")

    if (length(significantModules) == 0) {
        cat("No significant modules found for this trait.\n")
        return(NULL)
    }

    ## Create directory for this trait -----------------------------------------
    traitDirName <- gsub("[^[:alnum:]]", "_", traitName)
    traitDir <- file.path(outputDir, "module_analysis", traitDirName)
    dir.create(traitDir, showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(traitDir, "csv"), showWarnings = FALSE)
    dir.create(file.path(traitDir, "plots"), showWarnings = FALSE)

    ## Initialize list to store results for this trait -------------------------
    traitResults <- list()

    ## Loop through each significant module ------------------------------------
    for (moduleName in significantModules) {
        cat("\n--- Module:", moduleName, "for trait:", traitName, "---\n")

        ## Remove "ME" prefix for color name
        colorName <- gsub("ME", "", moduleName)

        ## Get genes in this module
        moduleGenes <- names(blockwiseNetwork$colors)[blockwiseNetwork$colors == colorName]
        cat("Genes in module:", length(moduleGenes), "\n")

        ## Calculate Gene Significance (GS) for this specific trait ------------
        traitValues <- traitMatrix[, traitColumnIndex]

        ## Check if genes exist in normalizedCounts
        availableGenes <- moduleGenes[moduleGenes %in% colnames(normalizedCounts)]
        if (length(availableGenes) < length(moduleGenes)) {
            cat("Warning:", length(moduleGenes) - length(availableGenes),
                "genes not found in expression matrix\n")
        }

        geneTraitCorrelation <- cor(
            normalizedCounts[, availableGenes],
            traitValues,
            use = "pairwise.complete.obs"
        )

        ## Calculate Module Membership (MM) ------------------------------------
        geneModuleMembership <- cor(
            normalizedCounts[, availableGenes],
            moduleEigengenes[, moduleName],
            use = "pairwise.complete.obs"
        )

        ## Combine results -----------------------------------------------------
        moduleData <- data.frame(
            gene_id = availableGenes,
            module_membership = as.numeric(geneModuleMembership),
            gene_significance = as.numeric(geneTraitCorrelation),
            stringsAsFactors = FALSE
        )

        ## Filter by user thresholds -------------------------------------------
        filteredGenes <- moduleData[
            abs(moduleData$gene_significance) > minGeneSignificance &
            abs(moduleData$module_membership) > minModuleMembership,
        ]

        cat("Genes passing thresholds:", nrow(filteredGenes), "\n")

        ## Initialize results for this module ----------------------------------
        moduleResults <- list(
            hub_genes = filteredGenes$gene_id,
            module_genes = availableGenes,
            top_hubgenes_expression_matrix = NULL
        )

        if (nrow(filteredGenes) >= 3) {
            ## Get top N hub genes by connectivity -----------------------------
            filteredGenesClean <- filteredGenes[complete.cases(filteredGenes), ]
            genesForAdjacency <- filteredGenesClean$gene_id

            if (length(genesForAdjacency) >= 3) {
                ## Check if all genes exist in normalizedCounts ----------------
                genesInMatrix <- genesForAdjacency[genesForAdjacency %in% colnames(normalizedCounts)]

                if (length(genesInMatrix) < length(genesForAdjacency)) {
                    cat("Warning:", length(genesForAdjacency) - length(genesInMatrix),
                        "genes not available for adjacency calculation\n")
                }

                if (length(genesInMatrix) >= 3) {
                    tryCatch({
                        adjacencyMatrix <- adjacency(
                            normalizedCounts[, genesInMatrix],
                            power = softPower,
                            type = "signed"
                        )

                        if (any(is.na(adjacencyMatrix))) {
                            cat("Warning: NA values in adjacency matrix, using na.rm=TRUE\n")
                            connectivityScores <- rowSums(adjacencyMatrix, na.rm = TRUE) - 1
                        } else {
                            connectivityScores <- rowSums(adjacencyMatrix) - 1
                        }

                        filteredGenesClean$connectivity <- connectivityScores[genesInMatrix]

                        ## Sort by connectivity and get top N -----------------
                        filteredGenesClean <- filteredGenesClean[
                            order(filteredGenesClean$connectivity, decreasing = TRUE),
                        ]

                        topHubGenes <- head(filteredGenesClean, nTopGenes)
                        moduleResults$top_hubgenes_expression_matrix <- t(
                            normalizedCounts[, topHubGenes$gene_id]
                        )

                        ## Save hub genes to CSV -------------------------------
                        write.csv(
                            moduleData,
                            file.path(traitDir, "csv", paste0(colorName, "_genes.csv")),
                            row.names = FALSE
                        )
                        write.csv(
                            filteredGenesClean,
                            file.path(traitDir, "csv", paste0(colorName, "_hub_genes.csv")),
                            row.names = FALSE
                        )
                        write.csv(
                            topHubGenes,
                            file.path(traitDir, "csv", paste0(colorName, "_top", nTopGenes, "_hub_genes.csv")),
                            row.names = FALSE
                        )

                        ## Plot MM vs GS for this module -----------------------
                        pdf(
                            file.path(traitDir, "plots", paste0("MM_vs_GS_", colorName, ".pdf")),
                            width = 8, height = 6
                        )
                        verboseScatterplot(
                            abs(moduleData$module_membership),
                            abs(moduleData$gene_significance),
                            xlab = paste("Module Membership in", colorName, "module"),
                            ylab = paste("Gene significance for", traitName),
                            main = paste("Module membership vs. gene significance\n"),
                            cex.main = 1.2,
                            cex.lab = 1.2,
                            cex.axis = 1.2,
                            col = colorName,
                            abline.color = "red",
                            abline = TRUE
                        )
                        dev.off()

                        cat("Saved hub genes CSV and MM vs GS plot\n")

                    }, error = function(error) {
                        cat("Error in adjacency calculation:", error$message, "\n")
                        cat("Skipping connectivity calculation for this module.\n")
                    })
                } else {
                    cat("Not enough valid genes for adjacency calculation (need at least 3).\n")
                }
            } else {
                cat("Not enough genes after removing NAs for adjacency calculation.\n")
            }
        } else {
            cat("Not enough genes for PPI analysis (need at least 3).\n")
        }

        ## Store module results in trait results -------------------------------
        traitResults[[colorName]] <- moduleResults
    }

    ## Return the results for this trait ---------------------------------------
    return(traitResults)
}

## =============================================================================
## Process each trait column
## =============================================================================

allTraitResults <- list()

## Get all trait columns
traitColumns <- colnames(traitMatrix)

## Process each trait column
for (i in seq_along(traitColumns)) {
    traitName <- traitColumns[i]
    result <- processSignificantModules(
        traitColumnIndex = i,
        traitName = traitName,
        nTopGenes = topNGenes
    )

    if (!is.null(result) && length(result) > 0) {
        allTraitResults[[traitName]] <- result
    }
}

## Save all results ------------------------------------------------------------
if (length(allTraitResults) > 0) {
    saveRDS(allTraitResults, file.path(outputDir, "all_traits_modules_hubgenes.rds"))
} else {
    saveRDS(NULL, file.path(outputDir, "all_traits_modules_hubgenes.rds"))
}
