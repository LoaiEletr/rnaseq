#!/usr/bin/env Rscript

## =============================================================================
## Script:  run_tximport.R
## Purpose: Import transcript quantification data using tximport
## Usage:   Rscript run_tximport.R <quant_dir1> <quant_dir2> ... <quant_type> <tx2gene_table>
## Example: Rscript run_tximport.R salmon_results1 salmon_results2 salmon tx2gene.tsv
## =============================================================================

## Load required packages ------------------------------------------------------
library(tximport)
library(rhdf5)
library(jsonlite)

## Parse command line arguments ------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop("Usage: run_tximport.R <quant_dir1> <quant_dir2> ... <quant_type> <tx2gene_table>",
         call. = FALSE)
}

## Extract arguments -----------------------------------------------------------
## All but last two arguments are quantification directories
quantDirs <- args[1:(length(args) - 2)]
## Second last argument is quantification type (salmon/kallisto)
quantType <- args[length(args) - 1]
## Last argument is tx2gene mapping file
tx2geneFile <- args[length(args)]

## =============================================================================
## Prepare file paths based on quantification type
## =============================================================================

## Determine file names based on quantification tool
if (quantType == "salmon") {
    abundanceFiles <- file.path(quantDirs, "quant.sf")
} else if (quantType == "kallisto") {
    abundanceFiles <- file.path(quantDirs, "abundance.tsv")
} else {
    stop("Unsupported quantification type: ", quantType,
         ". Supported types: salmon, kallisto", call. = FALSE)
}

## Verify all input files exist -----------------------------------------------
missingFiles <- abundanceFiles[!file.exists(abundanceFiles)]
if (length(missingFiles) > 0) {
    stop("Missing quantification files:\n", paste(missingFiles, collapse = "\n"),
         call. = FALSE)
}

## =============================================================================
## Read transcript-to-gene mapping
## =============================================================================

## Read tx2gene mapping file
tx2geneData <- read.delim(tx2geneFile, header = TRUE)

## Verify tx2gene has required columns
requiredCols <- c("target_id", "gene_name")
if (!all(requiredCols %in% colnames(tx2geneData))) {
    stop("tx2gene file must contain columns: ", paste(requiredCols, collapse = ", "),
         call. = FALSE)
}

## =============================================================================
## Add sample names to file vector
## =============================================================================

## Assign sample names as names of the abundanceFiles vector
sampleNames <- basename(quantDirs)
names(abundanceFiles) <- sampleNames

## =============================================================================
## Run tximport with named files
## =============================================================================

## Import quantification data (files have sample names as vector names)
txImportResult <- tximport(
    files = abundanceFiles,
    type = quantType,
    countsFromAbundance = "lengthScaledTPM",
    tx2gene = tx2geneData,
    txOut = FALSE,
    ignoreTxVersion = TRUE
)

## =============================================================================
## Save results
## =============================================================================

outputFile <- "tximport_results.rds"
saveRDS(txImportResult, file = outputFile)
