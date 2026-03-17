#!/usr/bin/env Rscript

## =============================================================================
## Script:  run_go_enrichment.R
## Purpose: Perform GO enrichment analysis for DE genes, WGCNA modules, MaSigPro clusters, DIU, and DEU genes
## Usage:   Rscript run_go_enrichment.R <degFile> <pvalThreshold> <lfcThreshold> <speciesName> <wgcnaList> <masigproList> <diuGenes> <deuGenes> <ntopProcesses>
## Example: Rscript run_go_enrichment.R deg.rds 0.05 1.0 human wgcna.rds masigpro.rds diu.rds deu.rds 20
## =============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(org.Sc.sgd.db)
library(org.Dm.eg.db)
library(org.Dr.eg.db)
library(org.Ce.eg.db)
library(org.At.tair.db)
library(org.Gg.eg.db)
library(org.Bt.eg.db)
library(org.Ss.eg.db)
library(org.Cf.eg.db)
library(org.Mmu.eg.db)
library(purrr)
library(biomaRt)
library(ggplot2)
library(tibble)
library(dplyr)

## =============================================================================
## Helper Functions
## =============================================================================

writeReasonPng <- function(file, reason) {
    png(file, width = 800, height = 600)
    plot.new()
    text(0.5, 0.5, reason, cex = 1.1)
    dev.off()
}

writeReasonCsv <- function(file, reason) {
    write.csv(
        data.frame(status = "NO_OUTPUT", reason = reason, time = Sys.time()),
        file, row.names = FALSE
    )
}

## =============================================================================
## Argument Parsing
## =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 9) {
    stop(
        "Usage:\n",
        "run_go_enrichment.R <degFile> <pvalThreshold> <lfcThreshold> <speciesName> <wgcnaList> <masigproList> <diuGenes> <deuGenes> <ntopProcesses>",
        call. = FALSE
    )
}

## Validate and set arguments
degFile          <- args[1]
pvalThreshold    <- as.numeric(args[2])
lfcThreshold     <- ifelse(args[3] == "null", 0, as.numeric(args[3]))
speciesName      <- args[4]
wgcnaListFile    <- args[5]
masigproListFile <- args[6]
diuGenesFile     <- args[7]
deuGenesFile     <- args[8]
ntopProcesses    <- as.numeric(args[9])

## =============================================================================
## Species Database
## =============================================================================

speciesDb <- list(
    human = list(orgdb = "org.Hs.eg.db", biomart = "hsapiens_gene_ensembl"),
    mouse = list(orgdb = "org.Mm.eg.db", biomart = "mmusculus_gene_ensembl"),
    rat   = list(orgdb = "org.Rn.eg.db", biomart = "rnorvegicus_gene_ensembl"),
    yeast = list(orgdb = "org.Sc.sgd.db", biomart = "scerevisiae_gene_ensembl"),
    fruitfly = list(orgdb = "org.Dm.eg.db", biomart = "dmelanogaster_gene_ensembl"),
    zebrafish = list(orgdb = "org.Dr.eg.db", biomart = "drerio_gene_ensembl"),
    worm = list(orgdb = "org.Ce.eg.db", biomart = "celegans_gene_ensembl"),
    arabidopsis = list(orgdb = "org.At.tair.db", biomart = "athaliana_eg_gene"),
    chicken = list(orgdb = "org.Gg.eg.db", biomart = "ggallus_gene_ensembl"),
    cow = list(orgdb = "org.Bt.eg.db", biomart = "btaurus_gene_ensembl"),
    pig = list(orgdb = "org.Ss.eg.db", biomart = "sscrofa_gene_ensembl"),
    dog = list(orgdb = "org.Cf.eg.db", biomart = "clfamiliaris_gene_ensembl"),
    monkey = list(orgdb = "org.Mmu.eg.db", biomart = "mmulatta_gene_ensembl")
)

if (!speciesName %in% names(speciesDb)) {
    stop("Unsupported species: ", speciesName)
}

speciesInfo <- speciesDb[[speciesName]]
dir.create("GO_results", showWarnings = FALSE)

## =============================================================================
## Load OrgDb for the Species
## =============================================================================

orgDbObject <- NULL
if (!is.null(speciesInfo$orgdb)) {
    tryCatch({
        ## Load the OrgDb package and get the object
        orgDbObject <- get(speciesInfo$orgdb)
        cat("Successfully loaded OrgDb:", speciesInfo$orgdb, "\n")

        ## Check available keytypes
        cat("Available keytypes in", speciesInfo$orgdb, ":\n")
        print(keytypes(orgDbObject))
    }, error = function(e) {
        cat("Warning: Could not load OrgDb", speciesInfo$orgdb, ":", e$message, "\n")
        orgDbObject <- NULL
    })
}

## =============================================================================
## Symbol → Entrez Conversion
## =============================================================================

convertSymbols <- function(genes) {

    ## Function to check if input contains Ensembl IDs
    isEnsembl <- function(geneVec) {
        sampleSize <- min(10, length(geneVec))
        sampleGenes <- geneVec[1:sampleSize]

        ## Check for Ensembl patterns
        ensemblPatterns <- c(
            "^ENS[A-Z]*G\\d+",  ## ENSG00000123456, ENSMUSG00000123456, etc.
            "^ENS[A-Z]*T\\d+",  ## ENST00000123456 (transcripts)
            "^[A-Z]{2}\\d{6,}"  ## Alternative Ensembl patterns like FBgn0001234
        )

        any(sapply(ensemblPatterns, function(pattern) {
            all(grepl(pattern, sampleGenes, ignore.case = TRUE))
        }))
    }

    ## Function to check if input contains Entrez IDs
    isEntrez <- function(geneVec) {
        sampleSize <- min(10, length(geneVec))
        sampleGenes <- geneVec[1:sampleSize]
        all(grepl("^\\d+$", sampleGenes))
    }

    ## Check if species is human, mouse, or rat - use OrgDb for these
    if (speciesName %in% c("human", "mouse", "rat")) {
        cat("Using OrgDb for", speciesName, "\n")

        if (!is.null(orgDbObject)) {
            result <- tryCatch({
                ## Check what type of IDs we have
                if (isEntrez(genes)) {
                    cat("Genes appear to be ENTREZ IDs\n")
                    return(genes)  ## Return as is
                } else if (isEnsembl(genes)) {
                    ## Convert ENSEMBL to ENTREZID
                    cat("Converting ENSEMBL to ENTREZID using OrgDb\n")

                    ## Determine the correct keytype for ENSEMBL IDs
                    keytypesList <- keytypes(orgDbObject)

                    ## Try different possible ENSEMBL keytypes
                    possibleEnsemblKeys <- c("ENSEMBL", "ENSEMBLTRANS", "ENSEMBLPROT", "ENSEMBLID")
                    ensemblKey <- NULL

                    for (key in possibleEnsemblKeys) {
                        if (key %in% keytypesList) {
                            ensemblKey <- key
                            break
                        }
                    }

                    if (is.null(ensemblKey)) {
                        cat("No ENSEMBL keytype found in OrgDb\n")
                        return(genes)
                    }

                    cat("Using keytype:", ensemblKey, "\n")

                    ## Map ENSEMBL to ENTREZID
                    conversion <- mapIds(orgDbObject,
                                       keys = genes,
                                       column = "ENTREZID",
                                       keytype = ensemblKey,
                                       multiVals = "first")

                    ## Remove NAs
                    entrezIds <- as.character(conversion[!is.na(conversion)])

                    if (length(entrezIds) > 0) {
                        return(entrezIds)
                    } else {
                        cat("No ENTREZ IDs found for ENSEMBL IDs\n")
                        return(genes)
                    }
                } else {
                    ## Assume SYMBOL and convert to ENTREZID
                    cat("Converting SYMBOL to ENTREZID using OrgDb\n")
                    conversion <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgDbObject)

                    ## Extract ENTREZID column as vector
                    if ("ENTREZID" %in% colnames(conversion)) {
                        return(conversion$ENTREZID)
                    } else {
                        cat("No ENTREZID column in conversion result\n")
                        return(genes)  ## Fallback: return input
                    }
                }
            }, error = function(e) {
                cat("OrgDb conversion failed:", e$message, "\n")
                return(genes)  ## Return input as fallback
            })
        }
    }

    ## For all other species, use biomaRt
    cat("Using biomaRt for", speciesName, "\n")

    ## Connect to BioMart
    mart <- if (grepl("_eg_", speciesInfo$biomart) || speciesName == "arabidopsis") {
        tryCatch({
            useMart("ensemblgenomes", dataset = speciesInfo$biomart)
        }, error = function(e) {
            tryCatch({
                useMart("plants_mart", dataset = speciesInfo$biomart, host = "https://may2024-plants.ensembl.org")
            }, error = function(e2) {
                useMart("ensemblgenomes_mart", dataset = speciesInfo$biomart, host = "https://ensemblgenomes.org")
            })
        })
    } else {
        useMart("ensembl", dataset = speciesInfo$biomart, host = "https://may2025.archive.ensembl.org")
    }

    if (!inherits(mart, "Mart")) {
        return(genes)
    }

    ## Check available attributes
    availableAttrs <- listAttributes(mart)$name

    ## Find ENTREZ ID attribute
    entrezAttr <- NULL
    possibleEntrezNames <- c("entrezgene_id", "entrezgene", "entrez_gene_id", "ncbi_gene_id")

    for (attrName in possibleEntrezNames) {
        if (attrName %in% availableAttrs) {
            entrezAttr <- attrName
            break
        }
    }

    ## If still not found, try pattern matching
    if (is.null(entrezAttr)) {
        entrezMatches <- grep("entrez|ncbi", availableAttrs, value = TRUE, ignore.case = TRUE)
        if (length(entrezMatches) > 0) {
            entrezAttr <- entrezMatches[1]
        }
    }

    ## If no ENTREZ attribute found, return input
    if (is.null(entrezAttr)) {
        cat("No ENTREZ ID attribute found for", speciesName, "\n")
        return(genes)
    }

    cat("Using attribute:", entrezAttr, "for ENTREZ IDs\n")

    ## Determine input type
    if (isEntrez(genes)) {
        cat("Input appears to be ENTREZ IDs\n")
        return(genes)
    } else if (isEnsembl(genes)) {
        cat("Input appears to be ENSEMBL IDs\n")

        ## Find Ensembl attribute
        ensemblAttr <- NULL
        possibleEnsemblNames <- c("ensembl_gene_id", "ensembl_transcript_id",
                                 "ensembl_peptide_id", "ensembl_gene_id_version",
                                 "ensembl_transcript_id_version")

        for (attrName in possibleEnsemblNames) {
            if (attrName %in% availableAttrs) {
                ensemblAttr <- attrName
                break
            }
        }

        ## Try pattern matching for Ensembl
        if (is.null(ensemblAttr)) {
            ensemblMatches <- grep("ensembl", availableAttrs, value = TRUE, ignore.case = TRUE)
            if (length(ensemblMatches) > 0) {
                ensemblAttr <- ensemblMatches[1]
            }
        }

        if (is.null(ensemblAttr)) {
            cat("No ENSEMBL attribute found\n")
            return(genes)
        }

        ## Convert ENSEMBL to ENTREZ IDs
        result <- tryCatch({
            chunkSize <- 500
            allEntrez <- c()

            for (i in seq(1, length(genes), chunkSize)) {
                chunk <- genes[i:min(i + chunkSize - 1, length(genes))]

                bmResult <- getBM(
                    attributes = c(entrezAttr),
                    filters = ensemblAttr,
                    values = chunk,
                    mart = mart
                )

                if (nrow(bmResult) > 0) {
                    allEntrez <- c(allEntrez, bmResult[[entrezAttr]])
                }
            }

            if (length(allEntrez) > 0) {
                return(allEntrez)
            } else {
                cat("No ENTREZ IDs found for ENSEMBL IDs\n")
                return(genes)
            }
        }, error = function(e) {
            cat("BioMart conversion failed:", e$message, "\n")
            return(genes)
        })

        return(result)

    } else {
        ## Assume input is SYMBOL
        cat("Input appears to be gene symbols\n")

        ## Find symbol attribute for conversion
        symbolAttr <- NULL
        possibleSymbolNames <- c("external_gene_name", "gene_name", "external_gene_id",
                               "hgnc_symbol", "mgi_symbol")

        for (attrName in possibleSymbolNames) {
            if (attrName %in% availableAttrs) {
                symbolAttr <- attrName
                break
            }
        }

        if (is.null(symbolAttr)) {
            cat("No symbol attribute found, cannot convert\n")
            return(genes)
        }

        ## Convert SYMBOL to ENTREZ IDs
        result <- tryCatch({
            chunkSize <- 500
            allEntrez <- c()

            for (i in seq(1, length(genes), chunkSize)) {
                chunk <- genes[i:min(i + chunkSize - 1, length(genes))]

                bmResult <- getBM(
                    attributes = c(entrezAttr),
                    filters = symbolAttr,
                    values = chunk,
                    mart = mart
                )

                if (nrow(bmResult) > 0) {
                    allEntrez <- c(allEntrez, bmResult[[entrezAttr]])
                }
            }

            if (length(allEntrez) > 0) {
                return(allEntrez)
            } else {
                cat("No ENTREZ IDs found for symbols\n")
                return(genes)
            }
        }, error = function(e) {
            cat("BioMart conversion failed:", e$message, "\n")
            return(genes)
        })

        return(result)
    }
}

## =============================================================================
## Run GO Enrichment
## =============================================================================

runGoEnrichment <- function(symbols, label) {

    ## Create directory structure based on label
    createDirsFromLabel <- function(label) {
        parts <- strsplit(label, "_")[[1]]

        if (length(parts) == 1) {
            dir.create(file.path("GO_results", label), showWarnings = FALSE, recursive = TRUE)
            baseDir <- file.path("GO_results", label)
            plotDir <- file.path(baseDir, "plots")
            csvDir <- file.path(baseDir, "csv")
        } else if (length(parts) == 2) {
            dir.create(file.path("GO_results", parts[1]), showWarnings = FALSE, recursive = TRUE)
            dir.create(file.path("GO_results", parts[1], parts[2]), showWarnings = FALSE, recursive = TRUE)
            baseDir <- file.path("GO_results", parts[1], parts[2])
            plotDir <- file.path(baseDir, "plots")
            csvDir <- file.path(baseDir, "csv")
        } else {
            dirPath <- file.path("GO_results", paste(parts[-length(parts)], collapse = "/"))
            subdir <- parts[length(parts)]
            dir.create(file.path(dirPath, subdir), showWarnings = FALSE, recursive = TRUE)
            baseDir <- file.path(dirPath, subdir)
            plotDir <- file.path(baseDir, "plots")
            csvDir <- file.path(baseDir, "csv")
        }

        dir.create(plotDir, showWarnings = FALSE, recursive = TRUE)
        dir.create(csvDir, showWarnings = FALSE, recursive = TRUE)

        return(list(baseDir = baseDir, plotDir = plotDir, csvDir = csvDir))
    }

    ## Get directory structure
    dirs <- createDirsFromLabel(label)
    plotDir <- dirs$plotDir
    csvDir <- dirs$csvDir

    if (is.null(symbols) || length(symbols) < 5) {
        writeReasonPng(
            file.path(plotDir, paste0(label, "_GO.png")),
            paste0("Not enough genes for GO enrichment: ", label)
        )
        return()
    }

    ## Check if OrgDb is available for enrichment
    if (is.null(orgDbObject)) {
        writeReasonPng(
            file.path(plotDir, paste0(label, "_GO.png")),
            paste0("No OrgDb available for ", speciesName, ". Cannot run GO enrichment.")
        )
        return()
    }

    for (ontology in c("BP", "CC", "MF")) {
        ## Use OrgDb for enrichmentGO
        goResult <- tryCatch({
            enrichGO(
                gene = symbols,
                OrgDb = orgDbObject,
                ont = ontology,
                pvalueCutoff = pvalThreshold
            )
        }, error = function(e) {
            cat("Error in enrichGO for", label, "(", ontology, "):", e$message, "\n")
            NULL
        })

        if (!is.null(goResult) && nrow(goResult) > 0) {
            tryCatch({
                ggsave(
                    file.path(plotDir, paste0(label, "_GO_", ontology, ".pdf")),
                    dotplot(goResult, showCategory = ntopProcesses) +
                        ggtitle(paste("GO", ontology, "Enrichment -", label)),
                    width = 12, height = 10
                )
                write.csv(goResult@result,
                        file.path(csvDir, paste0(label, "_GO_", ontology, ".csv")),
                        row.names = FALSE)
                cat("Successfully generated GO", ontology, "for", label, "\n")
            }, error = function(e) {
                cat("Error saving results for", label, "(", ontology, "):", e$message, "\n")
            })
        } else {
            writeReasonPng(
                file.path(plotDir, paste0(label, "_GO_", ontology, ".png")),
                paste0("No significant GO ", ontology, " terms for ", label)
            )
            cat("No significant GO", ontology, "terms for", label, "\n")
        }
    }
}

## =============================================================================
## Main Execution
## =============================================================================

## Basic DEG validation
if (degFile == "null") {
    writeReasonPng("GO_results/GO_DE.png", "No DE genes detected")
} else {
    degData <- readRDS(degFile)
    if (is.null(degData)) {
        writeReasonPng("GO_results/GO_DE.png", "No DE genes detected")
        quit(save = "no", status = 0)
    }

    upGenes   <- rownames(degData)[degData$logFC > lfcThreshold & degData$adj.P.Val < pvalThreshold]
    downGenes <- rownames(degData)[degData$logFC < -lfcThreshold & degData$adj.P.Val < pvalThreshold]

    upSymbols   <- convertSymbols(upGenes)
    downSymbols <- convertSymbols(downGenes)

    runGoEnrichment(upSymbols, "DEGs_Up")
    runGoEnrichment(downSymbols, "DEGs_Down")
}

## WGCNA modules
if (wgcnaListFile == "null") {
    writeReasonPng("GO_results/GO_WGCNA.png", "No genes detected from WGCNA modules")
} else {
    wgcnaData <- readRDS(wgcnaListFile)
    if (is.null(wgcnaData) || length(wgcnaData) == 0) {
        writeReasonPng("GO_results/GO_WGCNA.png", "No genes detected from WGCNA modules")
        quit(save = "no", status = 0)
    }

    for (traitName in names(wgcnaData)) {
        traitModules <- wgcnaData[[traitName]]

        if (is.null(traitModules) || length(traitModules) == 0) {
            next
        }

        ## Process each module for this trait
        for (moduleName in names(traitModules)) {
            moduleData <- traitModules[[moduleName]]

            ## Get the module genes
            if (is.null(moduleData$module_genes) || length(moduleData$module_genes) == 0) {
                next
            }

            moduleGenes <- moduleData$module_genes
            symbolsConverted <- convertSymbols(moduleGenes)
            label <- paste0(traitName, "_", moduleName)

            runGoEnrichment(symbolsConverted, label)
        }
    }
}

## MaSigPro clusters
if (masigproListFile == "null") {
    writeReasonPng("GO_results/GO_masigpro.png", "No genes detected from MaSigpro clusters")
} else {
    clusters <- readRDS(masigproListFile)
    if (is.null(clusters) || length(clusters) == 0) {
        writeReasonPng("GO_results/GO_masigpro.png", "No genes detected from MaSigpro clusters")
        quit(save = "no", status = 0)
    }

    for (clusterNum in sort(unique(clusters$cut))) {
        clusterGenes <- names(clusters$cut)[clusters$cut == clusterNum]

        if (is.null(clusterGenes) || length(clusterGenes) == 0) {
            next
        }

        symbolsConverted <- convertSymbols(clusterGenes)
        label <- paste0("Cluster", clusterNum)

        runGoEnrichment(symbolsConverted, label)
    }
}

## DIU genes
if (diuGenesFile == "null") {
    cat("No DIU genes file provided, skipping DIU GO enrichment\n")
} else {
    diuGenes <- readRDS(diuGenesFile)
    if (is.null(diuGenes) || length(diuGenes) == 0) {
        writeReasonPng("GO_results/GO_DIU.png", "No genes detected from DIU genes")
        cat("No DIU genes found\n")
    } else {
        symbolsConverted <- convertSymbols(diuGenes)
        runGoEnrichment(symbolsConverted, "DIU")
        cat("DIU GO enrichment completed\n")
    }
}

## DEU genes
if (deuGenesFile == "null") {
    cat("No DEU genes file provided, skipping DEU GO enrichment\n")
} else {
    deuGenes <- readRDS(deuGenesFile)
    if (is.null(deuGenes) || length(deuGenes) == 0) {
        writeReasonPng("GO_results/GO_DEU.png", "No genes detected from DEU genes")
        cat("No DEU genes found\n")
    } else {
        symbolsConverted <- convertSymbols(deuGenes)
        runGoEnrichment(symbolsConverted, "DEU")
        cat("DEU GO enrichment completed\n")
    }
}
