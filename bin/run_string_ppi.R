#!/usr/bin/env Rscript

## =============================================================================
## Script:  run_string_ppi.R
## Purpose: Protein-Protein Interaction analysis using STRING database
## Usage:   Rscript run_string_ppi.R <wgcna_results_file> <species_name> <score_threshold>
## Example: Rscript run_string_ppi.R wgcna_results.rds human 400
## Available species: human, mouse, rat, yeast, fruitfly, zebrafish, worm,
##                   arabidopsis, chicken, cow, pig, dog, monkey
## =============================================================================

## Load required packages ------------------------------------------------------
library(STRINGdb)

## =============================================================================
## Parse command line arguments
## =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop("Usage: run_string_ppi.R <wgcna_results_file> <species_name> <score_threshold>",
         call. = FALSE)
}

## Extract parameters with descriptive names -----------------------------------
wgcnaResultsFile <- args[1]
speciesName <- args[2]
scoreThreshold <- as.numeric(args[3])

## =============================================================================
## Define species database
## =============================================================================

speciesDatabase <- list(
    human = 9606,
    mouse = 10090,
    rat = 10116,
    yeast = 4932,
    fruitfly = 7227,
    zebrafish = 7955,
    worm = 6239,
    arabidopsis = 3702,
    chicken = 9031,
    cow = 9913,
    pig = 9823,
    dog = 9615,
    monkey = 9443
)

## Check if species is supported -----------------------------------------------
if (!speciesName %in% names(speciesDatabase)) {
    cat("Available species:", paste(names(speciesDatabase), collapse = ", "), "\n")
    stop("Unsupported species: ", speciesName, call. = FALSE)
}

speciesStringId <- speciesDatabase[[speciesName]]
cat("Using STRING database ID:", speciesStringId, "for species:", speciesName, "\n")

## =============================================================================
## Load WGCNA results
## =============================================================================

wgcnaResults <- readRDS(wgcnaResultsFile)

## Create output directory -----------------------------------------------------
outputDir <- "STRING_PPI"
dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)

## =============================================================================
## Process each trait and module from WGCNA results
## =============================================================================

## Check if wgcnaResults is valid ----------------------------------------------
if (is.null(wgcnaResults) || length(wgcnaResults) == 0) {
    cat("No significant modules found in WGCNA results.\n")
    quit(save = "no", status = 0)
}

## Initialize summary statistics -----------------------------------------------
summaryStatistics <- list(
    traits_processed = 0,
    modules_processed = 0,
    modules_with_ppi = 0,
    total_interactions = 0
)

## Loop through each trait in the results --------------------------------------
for (traitName in names(wgcnaResults)) {
    cat("\n=== Processing trait:", traitName, "===\n")
    summaryStatistics$traits_processed <- summaryStatistics$traits_processed + 1

    traitModules <- wgcnaResults[[traitName]]

    if (is.null(traitModules) || length(traitModules) == 0) {
        cat("No modules found for this trait.\n")
        next
    }

    ## Create trait directory --------------------------------------------------
    safeTraitName <- gsub("[^[:alnum:]_]", "_", traitName)
    traitDirectory <- file.path(outputDir, safeTraitName)
    dir.create(traitDirectory, showWarnings = FALSE, recursive = TRUE)

    ## Create subdirectories ---------------------------------------------------
    dir.create(
        file.path(traitDirectory, "STRING_cytoscape_edges"),
        showWarnings = FALSE,
        recursive = TRUE
    )
    dir.create(
        file.path(traitDirectory, "network_summaries"),
        showWarnings = FALSE,
        recursive = TRUE
    )

    ## Process each module for this trait --------------------------------------
    for (moduleName in names(traitModules)) {
        cat("\n--- Module:", moduleName, "---\n")
        summaryStatistics$modules_processed <- summaryStatistics$modules_processed + 1

        moduleData <- traitModules[[moduleName]]

        ## Check if module has hub genes ---------------------------------------
        if (is.null(moduleData$hub_genes) || length(moduleData$hub_genes) == 0) {
            cat("No hub genes available for PPI analysis.\n")
            next
        }

        genesToAnalyze <- unique(moduleData$hub_genes)  ## Remove duplicates
        cat("Number of hub genes to analyze:", length(genesToAnalyze), "\n")

        if (length(genesToAnalyze) < 2) {
            cat("Not enough genes for PPI analysis (need at least 2).\n")
            next
        }

        ## Perform STRING PPI analysis -----------------------------------------
        tryCatch({
            stringDatabase <- STRINGdb$new(
                species = speciesStringId,
                score_threshold = scoreThreshold
            )

            ## Map genes to STRING IDs -----------------------------------------
            hubGenesMapped <- stringDatabase$map(
                data.frame(gene = genesToAnalyze, stringsAsFactors = FALSE),
                "gene",
                removeUnmappedRows = TRUE
            )

            cat("Successfully mapped", nrow(hubGenesMapped), "genes to STRING database\n")

            if (nrow(hubGenesMapped) > 1) {
                ## Get interactions ---------------------------------------------
                ppiInteractions <- stringDatabase$get_interactions(hubGenesMapped$STRING_id)

                if (nrow(ppiInteractions) > 0) {
                    ## Remove duplicate interactions ----------------------------
                    ppiInteractions <- ppiInteractions[
                        !duplicated(ppiInteractions[, c("from", "to")]),
                    ]
                    cat("Found", nrow(ppiInteractions), "unique interactions\n")

                    ## Filter by score threshold --------------------------------
                    ppiInteractions <- ppiInteractions[
                        ppiInteractions$combined_score >= scoreThreshold,
                    ]
                    cat("After filtering (score >=", scoreThreshold, "):",
                        nrow(ppiInteractions), "interactions\n")

                    if (nrow(ppiInteractions) > 0) {
                        summaryStatistics$modules_with_ppi <- summaryStatistics$modules_with_ppi + 1
                        summaryStatistics$total_interactions <- summaryStatistics$total_interactions +
                            nrow(ppiInteractions)

                        ## Get protein information ------------------------------
                        proteinInformation <- stringDatabase$add_proteins_description(hubGenesMapped)

                        ## Create mappings --------------------------------------
                        proteinMapping <- setNames(
                            proteinInformation$preferred_name,
                            proteinInformation$STRING_id
                        )
                        stringToGene <- setNames(
                            hubGenesMapped$gene,
                            hubGenesMapped$STRING_id
                        )

                        ## Annotate interactions --------------------------------
                        annotatedInteractions <- data.frame(
                            from_gene = stringToGene[as.character(ppiInteractions$from)],
                            to_gene = stringToGene[as.character(ppiInteractions$to)],
                            combined_score = ppiInteractions$combined_score,
                            stringsAsFactors = FALSE
                        )

                        ## Remove NA rows ---------------------------------------
                        annotatedInteractions <- na.omit(annotatedInteractions)

                        ## Create edge list for Cytoscape -----------------------
                        edgeList <- data.frame(
                            from = proteinMapping[as.character(ppiInteractions$from)],
                            to = proteinMapping[as.character(ppiInteractions$to)],
                            combined_score = ppiInteractions$combined_score,
                            alt_from = annotatedInteractions$from_gene,
                            alt_to = annotatedInteractions$to_gene,
                            stringsAsFactors = FALSE
                        )

                        edgeList <- na.omit(edgeList)

                        if (nrow(edgeList) > 0) {
                            ## Save edge list for Cytoscape ---------------------
                            edgeBaseName <- paste0(moduleName, "_", safeTraitName, "_ppi")
                            edgeFilePath <- file.path(
                                traitDirectory,
                                "STRING_cytoscape_edges",
                                paste0(edgeBaseName, "_edges.csv")
                            )

                            write.csv(edgeList, edgeFilePath, row.names = FALSE)

                            ## Create network summary ---------------------------
                            networkSummary <- data.frame(
                                trait = traitName,
                                module = moduleName,
                                total_genes = length(genesToAnalyze),
                                mapped_genes = nrow(hubGenesMapped),
                                total_interactions = nrow(edgeList),
                                min_score = min(edgeList$combined_score),
                                max_score = max(edgeList$combined_score),
                                mean_score = mean(edgeList$combined_score),
                                edge_files = paste0(basename(edgeFilePath)),
                                stringsAsFactors = FALSE
                            )

                            summaryFilePath <- file.path(
                                traitDirectory,
                                "network_summaries",
                                paste0(moduleName, "_network_summary.csv")
                            )
                            write.csv(networkSummary, summaryFilePath, row.names = FALSE)
                            cat("Network summary saved to:", summaryFilePath, "\n")

                        } else {
                            cat("No valid interactions after filtering.\n")
                        }
                    } else {
                        cat("No interactions found above score threshold.\n")
                    }
                } else {
                    cat("No PPI interactions found for genes.\n")
                }
            } else {
                cat("Not enough genes mapped to STRING database (need at least 2).\n")
            }

        }, error = function(error) {
            cat("ERROR in PPI analysis for module", moduleName, ":", error$message, "\n")
        })
    }
}

## =============================================================================
## Save summary statistics
## =============================================================================

summaryDataFrame <- data.frame(
    metric = names(summaryStatistics),
    value = unlist(summaryStatistics)
)

write.csv(
    summaryDataFrame,
    file.path(outputDir, "analysis_summary.csv"),
    row.names = FALSE
)
