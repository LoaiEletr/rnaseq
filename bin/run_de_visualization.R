#!/usr/bin/env Rscript

## =============================================================================
## Script:  run_de_visualization.R
## Purpose: Generate heatmaps (DE genes, MaSigPro clusters, WGCNA modules) and volcano plot
## Usage:   Rscript run_de_visualization.R <degFile> <samplesheetFile> <clusterFile> <pvalThreshold> <lfcThreshold> <vstFile>
## Example: Rscript run_de_visualization.R deg.rds samples.csv clusters.rds 0.05 1.0 vst.rds
## =============================================================================

library(pheatmap)
library(EnhancedVolcano)
library(RColorBrewer)

## =============================================================================
## Helper Functions (All camelCase)
## =============================================================================

## Create cluster heatmap
createClusterHeatmap <- function(geneList,
                                 vstMat,
                                 metaData,
                                 clusterId,
                                 outputPath,
                                 width = 10,
                                 height = 10,
                                 showRowNames = FALSE,
                                 showColNames = FALSE,
                                 clusterRows = TRUE,
                                 clusterCols = FALSE,
                                 fontSizeRow = 7,
                                 winsorizeLimit = 2) {

    genesFound <- intersect(geneList, rownames(vstMat))

    if (length(genesFound) < 2) {
        warning("Insufficient genes found for heatmap: ", length(genesFound))
        return(NULL)
    }

    ## Check if 'timepoint' column exists in metadata
    hasTimepoint <- "timepoint" %in% colnames(metaData)

    ## Order samples based on presence of timepoint
    if (hasTimepoint) {
        ## Order samples: Disease by timepoint, then Control by timepoint
        diseaseIdx <- which(metaData$condition != "control")
        controlIdx <- which(metaData$condition == "control")

        metaDisease <- metaData[diseaseIdx, , drop = FALSE]
        metaControl <- metaData[controlIdx, , drop = FALSE]

        ## Order disease samples by timepoint
        diseaseOrder <- order(metaDisease$timepoint)
        metaDisease <- metaDisease[diseaseOrder, , drop = FALSE]

        ## Order control samples by timepoint
        controlOrder <- order(metaControl$timepoint)
        metaControl <- metaControl[controlOrder, , drop = FALSE]

        metaOrdered <- rbind(metaDisease, metaControl)
    } else {
        ## If no timepoint, order samples: Disease first, then Control
        diseaseIdx <- which(metaData$condition != "control")
        controlIdx <- which(metaData$condition == "control")

        metaDisease <- metaData[diseaseIdx, , drop = FALSE]
        metaControl <- metaData[controlIdx, , drop = FALSE]

        metaOrdered <- rbind(metaDisease, metaControl)
    }

    sampleOrder <- metaOrdered$sample_id
    vstSubset <- vstMat[genesFound, sampleOrder, drop = FALSE]

    ## Scale the data (z-score normalization) and handle NA/NaN/Inf
    vstZscore <- t(scale(t(vstSubset)))

    ## Replace infinite values with winsorized limits
    vstZscore[is.infinite(vstZscore) & vstZscore > 0] <- winsorizeLimit
    vstZscore[is.infinite(vstZscore) & vstZscore < 0] <- -winsorizeLimit

    ## Replace NA/NaN values with 0
    vstZscore[is.na(vstZscore) | is.nan(vstZscore)] <- 0

    ## Apply winsorization
    vstZscore[vstZscore > winsorizeLimit] <- winsorizeLimit
    vstZscore[vstZscore < -winsorizeLimit] <- -winsorizeLimit

    ## Check if there are still problematic values
    if (any(is.na(vstZscore)) || any(is.nan(vstZscore)) || any(is.infinite(vstZscore))) {
        warning("Matrix contains NA/NaN/Inf values after cleaning. Clustering may fail.")
    }

    ## Prepare annotation data frame based on available columns
    if (hasTimepoint) {
        annotationDf <- data.frame(
            Condition = metaOrdered$condition,
            Timepoint = factor(
                paste0("T", metaOrdered$timepoint),
                levels = unique(paste0("T", metaData$timepoint))
            ),
            row.names = metaOrdered$sample_id
        )
        colnames(annotationDf) <- c("Condition", "Timepoint")
    } else {
        annotationDf <- data.frame(
            Condition = metaOrdered$condition,
            row.names = metaOrdered$sample_id
        )
        colnames(annotationDf) <- "Condition"
    }

    ## Set up colors
    conditions <- unique(metaData$condition)
    nonControlConditions <- setdiff(conditions, "control")

    ## Create condition colors
    conditionColors <- c(control = "#1F78B4")
    if (length(nonControlConditions) > 0) {
        conditionColors <- c(conditionColors, "#E31A1C")
        names(conditionColors)[2] <- nonControlConditions[1]
    }

    annotationColors <- list(Condition = conditionColors)

    if (hasTimepoint) {
        timepointLevels <- unique(paste0("T", metaData$timepoint))
        timepointColors <- RColorBrewer::brewer.pal(
            length(unique(metaData$timepoint)),
            "Set2"
        )
        names(timepointColors) <- timepointLevels
        annotationColors[["Timepoint"]] <- timepointColors
    }

    heatmapColors <- grDevices::colorRampPalette(c(
        "#053061", "#2166AC", "#4393C3", "#92C5DE",
        "#D1E5F0", "#FDDBC7", "#F4A582", "#D6604D",
        "#B2182B", "#67001F"
    ))(100)

    quantileBreaks <- unique(quantile(as.vector(vstZscore),
                                     probs = seq(0, 1, length.out = 101),
                                     na.rm = TRUE))

    ## Gap position between Disease and Control groups (if applicable)
    if (hasTimepoint) {
        nDisease <- sum(metaOrdered$condition != "control")
        if (nDisease > 0) {
            gapPositions <- nDisease
        } else {
            gapPositions <- NULL
        }
    } else {
        gapPositions <- NULL
    }

    ## Add a safety check for clustering
    if (clusterRows) {
        ## Check if rows have variance
        rowVariances <- apply(vstZscore, 1, var, na.rm = TRUE)
        if (any(rowVariances == 0, na.rm = TRUE)) {
            warning("Some rows have zero variance. Disabling row clustering for these rows.")
        }
    }

    if (clusterCols) {
        ## Check if columns have variance
        colVariances <- apply(vstZscore, 2, var, na.rm = TRUE)
        if (any(colVariances == 0, na.rm = TRUE)) {
            warning("Some columns have zero variance. Disabling column clustering for these columns.")
        }
    }

    p <- pheatmap::pheatmap(
        vstZscore,
        color = heatmapColors,
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "euclidean",
        clustering_method = "complete",
        cluster_rows = clusterRows,
        cluster_cols = clusterCols,
        annotation_col = annotationDf,
        annotation_colors = annotationColors,
        show_rownames = showRowNames,
        show_colnames = showColNames,
        fontsize = 10,
        fontsize_row = fontSizeRow,
        fontsize_col = 8,
        border_color = NA,
        scale = "none",
        annotation_legend = TRUE,
        annotation_names_col = FALSE,
        legend = TRUE,
        breaks = quantileBreaks,
        gaps_col = gapPositions,
        main = paste0(clusterId, " (n=", length(geneList), " genes)"),
        silent = FALSE
    )

    pdf(outputPath, width = width, height = height)
    print(p)
    dev.off()

    return(p)
}

## =============================================================================
## Main Execution
## =============================================================================

## Parse command line arguments ------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
    stop("Usage: run_de_visualization.R <degFile> <samplesheetFile> <clusterFile> <pvalThreshold> <lfcThreshold> <vstFile>",
         call. = FALSE)
}

degFile <- args[1]
samplesheetFile <- args[2]
clusterFile <- args[3]
pvalThreshold <- as.numeric(args[4])
lfcThreshold <- as.numeric(args[5])
vstFile <- args[6]

## Load data -------------------------------------------------------------------
cat("Loading data...\n")

## Load DEG data for volcano plot
degData <- if (degFile != "null") readRDS(degFile) else NULL

## Load samplesheet
samplesheet <- read.csv(samplesheetFile)
if (is.null(samplesheet$timepoint)) {
    samplesheet <- samplesheet[, c("sample_id", "condition"), drop = FALSE]
} else {
    samplesheet <- samplesheet[, c("sample_id", "condition", "timepoint"), drop = FALSE]
}
rownames(samplesheet) <- samplesheet$sample_id

## Load VST data for heatmaps
vstData <- readRDS(vstFile)

## =============================================================================
## Generate Volcano Plot (for regular DE)
## =============================================================================

cat("Generating volcano plot for DE results...\n")
outputDir <- "de_visualization"
dir.create(outputDir, showWarnings = FALSE)

pdf(file.path(outputDir, "volcano_sig_genes.pdf"), width = 10, height = 8)

## CASE 1: NULL object
if (is.null(degData)) {
    plot.new()
    text(
        0.5, 0.5,
        "Volcano plot not generated:\nDifferential expression results are NULL.",
        cex = 1.3
    )

## CASE 2: Empty or invalid table
} else if (!is.data.frame(degData) ||
           nrow(degData) == 0 ||
           !all(c("logFC", "adj.P.Val") %in% colnames(degData))) {

    plot.new()
    text(
        0.5, 0.5,
        paste0(
            "Volcano plot not generated:\n",
            "Invalid or empty DE results table.\n\n",
            "Rows: ", ifelse(is.null(nrow(degData)), "NA", nrow(degData)), "\n",
            "Columns found: ", paste(colnames(degData), collapse = ", ")
        ),
        cex = 1.1
    )

## CASE 3: Normal volcano plot
} else {
    EnhancedVolcano(
        degData,
        lab = rownames(degData),
        x = "logFC",
        y = "adj.P.Val",
        xlab = bquote(~Log[2]~ "fold change"),
        ylab = bquote(~-Log[10]~ "adjusted p-value"),
        pCutoff = pvalThreshold,
        FCcutoff = lfcThreshold,
        pointSize = 3.0,
        labSize = 4.0,
        colAlpha = 0.8,
        legendPosition = "right",
        legendLabSize = 12,
        legendIconSize = 4.0,
        drawConnectors = TRUE,
        widthConnectors = 0.75,
        colConnectors = "black",
        max.overlaps = 15,
        caption = paste0("produced on ", Sys.time()),
        title = "Volcano Plot – All DE Genes"
    )
}

dev.off()
cat("✓ Volcano plot generated: de_visualization/volcano_sig_genes.pdf\n")

## =============================================================================
## Generate Heatmaps based on input type
## =============================================================================

cat("Generating heatmaps...\n")
heatmapDir <- file.path(outputDir, "heatmap")
dir.create(heatmapDir, showWarnings = FALSE)
errorPlotPath <- file.path(heatmapDir, "heatmap_failed.pdf")

## Check if VST data is WGCNA result (list structure)
if (is.list(vstData)) {

    cat("Processing WGCNA results...\n")
    wgcnaList <- vstData

    for (traitName in names(wgcnaList)) {

        traitModules <- wgcnaList[[traitName]]

        if (is.null(traitModules) || length(traitModules) == 0) {
            next
        }

        ## Process each module for this trait
        for (moduleName in names(traitModules)) {

            moduleData <- traitModules[[moduleName]]
            path <- file.path(heatmapDir, traitName)
            dir.create(path, recursive = TRUE, showWarnings = FALSE)

            ## Get the hub genes (note: original data uses snake_case)
            if (is.null(moduleData$hub_genes) || length(moduleData$hub_genes) == 0) {
                next
            }

            expressionMatrix <- moduleData$top_hubgenes_expression_matrix
            genes <- rownames(expressionMatrix)

            nGenes <- length(genes)
            showNames <- nGenes <= 50
            heightVal <- min(8, max(4, nGenes * 0.08))

            createClusterHeatmap(
                geneList = genes,
                vstMat = expressionMatrix,
                metaData = samplesheet,
                clusterId = paste0("hub genes for Module ", moduleName, " in ", traitName, " trait"),
                outputPath = file.path(path, paste0("modules_hubgenes_", moduleName, "_heatmap.pdf")),
                width = 10,
                height = heightVal,
                showRowNames = showNames,
                showColNames = FALSE,
                clusterRows = TRUE,
                clusterCols = FALSE,
                fontSizeRow = 7,
                winsorizeLimit = 3
            )
        }
    }
    cat("✓ WGCNA heatmaps generated\n")

} else if (clusterFile != "null" && file.exists(clusterFile)) {

    cat("Processing MaSigPro clusters...\n")
    expressionMatrix <- vstData
    clusters <- readRDS(clusterFile)
    masigproPath <- file.path(heatmapDir, "masigpro")
    dir.create(masigproPath, recursive = TRUE)

    for (cl in sort(unique(clusters$cut))) {
        clusterGenes <- names(clusters$cut)[clusters$cut == cl]

        nGenes <- length(clusterGenes)
        showNames <- nGenes <= 50
        heightVal <- min(8, max(4, nGenes * 0.08))

        createClusterHeatmap(
            geneList = clusterGenes,
            vstMat = expressionMatrix,
            metaData = samplesheet,
            clusterId = paste0("Cluster : ", cl),
            outputPath = file.path(masigproPath, paste0("cluster_", cl, "_heatmap.pdf")),
            width = 10,
            height = heightVal,
            showRowNames = showNames,
            showColNames = FALSE,
            clusterRows = TRUE,
            clusterCols = FALSE,
            fontSizeRow = 7,
            winsorizeLimit = 3
        )
    }
    cat("✓ MaSigPro heatmaps generated\n")

} else if (is.null(vstData) || is.vector(vstData) || nrow(vstData) < 2 || ncol(vstData) < 2) {

    pdf(errorPlotPath, width = 6, height = 4)
    plot.new()
    text(
        0.5, 0.5,
        paste0(
            "Heatmap not generated:\n",
            "Not enough data for clustering.\n",
            "Rows: ", ifelse(is.null(nrow(vstData)), "NA", nrow(vstData)),
            " | Columns: ", ifelse(is.null(ncol(vstData)), "NA", ncol(vstData))
        ),
        cex = 1.2
    )
    dev.off()
    cat("✗ Heatmap not generated: insufficient data\n")

} else {

    cat("Processing regular DE genes...\n")
    expressionMatrix <- vstData
    genes <- rownames(expressionMatrix)
    nGenes <- length(genes)
    showNames <- nGenes <= 50
    heightVal <- min(8, max(4, nGenes * 0.08))

    createClusterHeatmap(
        geneList = genes,
        vstMat = expressionMatrix,
        metaData = samplesheet,
        clusterId = paste0("Differentially Expressed Genes"),
        outputPath = file.path(heatmapDir, "DEGs_heatmap.pdf"),
        width = 10,
        height = heightVal,
        showRowNames = showNames,
        showColNames = FALSE,
        clusterRows = TRUE,
        clusterCols = FALSE,
        fontSizeRow = 7,
        winsorizeLimit = 3
    )
    cat("✓ DE genes heatmap generated\n")
}
