#!/usr/bin/env Rscript

## =============================================================================
## Script:  generate_tx2gene.R
## Purpose: Generate transcript-to-gene mapping file (tx2gene) for various species
## Usage:   Rscript generate_tx2gene.R <species_name>
## Supported species: human, mouse, rat, yeast, fruitfly, zebrafish, worm,
##                   arabidopsis, chicken, cow, pig, dog, monkey
## =============================================================================

## Load required packages ------------------------------------------------------
library(biomaRt)

## Parse command line arguments ------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
    stop("Usage: generate_tx2gene.R <species_name>", call. = FALSE)
}

speciesName <- tolower(args[1])  ## Convert to lowercase for case-insensitive matching

## =============================================================================
## BioMart dataset lookup for all species
## =============================================================================

## BioMart dataset mapping for all supported species
biomartLookup <- list(
    human = "hsapiens_gene_ensembl",
    mouse = "mmusculus_gene_ensembl",
    rat = "rnorvegicus_gene_ensembl",
    yeast = "scerevisiae_gene_ensembl",
    fruitfly = "dmelanogaster_gene_ensembl",
    zebrafish = "drerio_gene_ensembl",
    worm = "celegans_gene_ensembl",
    arabidopsis = "athaliana_eg_gene",
    chicken = "ggallus_gene_ensembl",
    cow = "btaurus_gene_ensembl",
    pig = "sscrofa_gene_ensembl",
    dog = "clfamiliaris_gene_ensembl",
    monkey = "mmulatta_gene_ensembl"
)

## Check if species is supported
if (!speciesName %in% names(biomartLookup)) {
    stop("Species '", speciesName, "' is not supported.", call. = FALSE)
}

dataset <- biomartLookup[[speciesName]]

## =============================================================================
## Connect to appropriate BioMart
## =============================================================================

mart <- tryCatch({
    if (grepl("_eg_", dataset) || speciesName == "arabidopsis") {
        ## Try ensemblgenomes for plants and non-vertebrates
        tryCatch({
            useMart("ensemblgenomes", dataset = dataset)
        }, error = function(e) {
            tryCatch({
                useMart("plants_mart", dataset = dataset, host = "https://plants.ensembl.org")
            }, error = function(e2) {
                useMart("ensemblgenomes_mart", dataset = dataset, host = "https://ensemblgenomes.org")
            })
        })
    } else {
        ## Standard Ensembl for vertebrates
        useMart("ensembl", dataset = dataset)
    }
}, error = function(e) {
    cat("Failed to connect to BioMart:", e$message, "\n")

    ## Try alternative Arabidopsis connection
    if (speciesName == "arabidopsis") {
        cat("Trying alternative Arabidopsis connection...\n")
        tryCatch({
            useMart("plants_mart", host = "https://plants.ensembl.org")
        }, error = function(e2) {
            stop("Could not connect to BioMart for ", speciesName, call. = FALSE)
        })
    } else {
        stop("Could not connect to BioMart for ", speciesName, call. = FALSE)
    }
})

## =============================================================================
## Retrieve data from BioMart
## =============================================================================

tx2geneData <- getBM(
    attributes = c("ensembl_transcript_id", "ensembl_gene_id"),
    mart = mart
)

## Rename columns using base R
colnames(tx2geneData) <- c("target_id", "gene_name")

## =============================================================================
## Data validation and output
## =============================================================================

## Check if data was retrieved successfully
if (nrow(tx2geneData) == 0) {
    stop("Failed to retrieve transcript-to-gene mapping for ", speciesName, call. = FALSE)
}

## Remove rows with missing gene names (base R)
tx2geneData <- tx2geneData[!is.na(tx2geneData$gene_name) & tx2geneData$gene_name != "", ]

## Write output file -----------------------------------------------------------
outputFile <- "tx2gene.tsv"
write.table(
    tx2geneData,
    file = outputFile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
