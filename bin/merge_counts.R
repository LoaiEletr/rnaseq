#!/usr/bin/env Rscript

## =============================================================================
## Script:  merge_counts.R
## Purpose: Merge multiple featureCounts output files into a single count matrix
## Usage:   Rscript merge_counts.R <count_file1> <count_file2> ... <count_fileN>
## =============================================================================

## Load required packages ------------------------------------------------------
library(tools)  # for file_path_sans_ext

## Parse command line arguments ------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Usage: merge_counts.R <count_file1> <count_file2> ... <count_fileN>",
         call. = FALSE)
}

countFiles <- args

## Helper function: clean sample name (remove everything after first dot) ------
cleanSampleName <- function(filename) {
    ## Remove file extension first
    nameNoExt <- file_path_sans_ext(filename)
    ## Remove everything after first dot (including the dot)
    cleanName <- sub("\\..*", "", nameNoExt)
    return(cleanName)
}

## Helper function: read and process a single featureCounts file ---------------
readFeatureCounts <- function(filePath) {
    ## Read the featureCounts output file
    ## Skip lines starting with '#' (header/comments)
    lines <- readLines(filePath)
    commentLines <- grep("^#", lines)
    if (length(commentLines) > 0) {
        lines <- lines[-commentLines]
    }

    ## Parse as data frame
    con <- textConnection(lines)
    countData <- read.delim(con, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    close(con)

    ## Extract and clean sample name from file path
    fileName <- basename(filePath)
    sampleName <- cleanSampleName(fileName)

    ## Keep only Geneid and counts (7th column)
    if (ncol(countData) >= 7) {
        result <- data.frame(
            Geneid = countData$Geneid,
            Counts = countData[, 7],
            stringsAsFactors = FALSE
        )
        colnames(result)[2] <- sampleName
        return(result)
    } else {
        stop("File ", filePath, " does not have enough columns")
    }
}

## Helper function: merge data frames by Geneid --------------------------------
mergeDataFrames <- function(dfList) {
    if (length(dfList) == 0) return(data.frame())
    if (length(dfList) == 1) return(dfList[[1]])

    ## Start with first data frame
    merged <- dfList[[1]]

    ## Merge with remaining data frames
    for (i in 2:length(dfList)) {
        merged <- merge(merged, dfList[[i]], by = "Geneid", all = TRUE)
    }

    ## Replace NA with 0 and sum duplicates
    merged[is.na(merged)] <- 0

    ## Aggregate duplicate Geneids (sum columns)
    geneIds <- merged$Geneid
    numericCols <- sapply(merged, is.numeric)

    ## Use aggregate to sum by Geneid
    result <- aggregate(merged[, numericCols, drop = FALSE],
                       by = list(Geneid = geneIds),
                       FUN = sum)

    return(result)
}

## Main processing pipeline ----------------------------------------------------
## 1. Read all count files
rawCounts <- lapply(countFiles, readFeatureCounts)

## 2. Merge all data frames by Geneid
mergedCounts <- mergeDataFrames(rawCounts)

## 3. Convert to count matrix (genes as rows, samples as columns)
if (nrow(mergedCounts) > 0) {
    countMatrix <- as.matrix(mergedCounts[, -1, drop = FALSE])
    rownames(countMatrix) <- mergedCounts$Geneid

    ## Clean column names (remove everything after dot)
    colnames(countMatrix) <- sapply(colnames(countMatrix), function(x) {
        if (grepl("\\.", x)) {
            return(sub("\\..*", "", x))
        } else {
            return(x)
        }
    })
} else {
    countMatrix <- matrix(nrow = 0, ncol = 0)
}

## Save output -----------------------------------------------------------------
outputFile <- "merged_counts.rds"
saveRDS(countMatrix, file = outputFile)
