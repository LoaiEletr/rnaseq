#!/usr/bin/env Rscript

## =============================================================================
## Script:  run_msigdb_gsea.R
## Purpose: Perform Gene Set Enrichment Analysis (GSEA) using MSigDB collections
## Usage:   Rscript run_msigdb_gsea.R <degFile> <pvalThreshold> <speciesName> <msigdbCategories> <nesThreshold> <padjGsea> <ntopProcesses> <rankMethod>
## Example: Rscript run_msigdb_gsea.R deg.rds 0.05 human C2,C5 1.5 0.1 20 t_stat
## =============================================================================

library(clusterProfiler)
library(enrichplot)
library(msigdbr)
library(dplyr)
library(tibble)
library(ggplot2)

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

if (length(args) < 8) {
    stop(
        "Usage:\n",
        "run_msigdb_gsea.R <degFile> <pvalThreshold> <speciesName> <msigdbCategories> <nesThreshold> <padjGsea> <ntopProcesses> <rankMethod>",
        call. = FALSE
    )
}

## =============================================================================
## Validate and Set Arguments
## =============================================================================

degFile          <- args[1]
pvalThreshold    <- as.numeric(args[2])
speciesName      <- args[3]
msigCategories   <- strsplit(args[4], ",")[[1]]
nesThreshold     <- as.numeric(args[5])
padjGsea         <- as.numeric(args[6])
ntopProcesses    <- as.numeric(args[7])
rankMethod       <- args[8]

## =============================================================================
## Species Database
## =============================================================================

speciesDb <- list(
    human = "Homo sapiens",
    mouse = "Mus musculus",
    rat   = "Rattus norvegicus",
    yeast = "Saccharomyces cerevisiae",
    fruitfly = "Drosophila melanogaster",
    zebrafish = "Danio rerio",
    worm = "Caenorhabditis elegans",
    arabidopsis = "Arabidopsis thaliana",
    chicken = "Gallus gallus",
    cow = "Bos taurus",
    pig = "Sus scrofa",
    dog = "Canis familiaris",
    monkey = "Macaca mulatta"
)

if (!speciesName %in% names(speciesDb)) {
    stop("Unsupported species: ", speciesName)
}

speciesLatin <- speciesDb[[speciesName]]

## =============================================================================
## Load DEG Data
## =============================================================================

degData <- readRDS(degFile)

## =============================================================================
## Create Output Directories
## =============================================================================

dir.create("GSEA_results/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("GSEA_results/csv", recursive = TRUE, showWarnings = FALSE)

## =============================================================================
## Prepare Ranked Gene List
## =============================================================================

prepareRankedGeneList <- function(degData, rankMethod) {
    if (rankMethod == "t_stat") {
        ranked <- degData$t
    } else if (rankMethod == "logfc") {
        ranked <- degData$logFC
    } else if (rankMethod == "signed_significance") {
        ranked <- degData$logFC * -log10(degData$P.Value)
    } else {
        stop("Invalid rank method. Choose from: t_stat, logfc, signed_significance")
    }

    names(ranked) <- rownames(degData)
    ranked <- sort(ranked, decreasing = TRUE)
    return(ranked)
}

rankedGenes <- prepareRankedGeneList(degData, rankMethod)

## =============================================================================
## Main GSEA Analysis
## =============================================================================

for (category in msigCategories) {

    ## Download MSigDB data for the species and category
    msigData <- msigdbr(species = speciesLatin, collection = category)

    if (nrow(msigData) == 0) {
        writeReasonPng(
            file.path("GSEA_results/plots", paste0("GSEA_", category, ".png")),
            paste0("No MSigDB genes for category: ", category)
        )
        next
    }

    ## Determine appropriate column for gene identifiers
    if (all(grepl("^ENS", names(rankedGenes)))) {
        termToGene <- msigData[, c("gs_name", "ensembl_gene")]
    } else {
        termToGene <- msigData[, c("gs_name", "gene_symbol")]
    }

    ## Run GSEA
    gseaResult <- GSEA(rankedGenes, TERM2GENE = termToGene, verbose = FALSE)

    if (is.null(gseaResult) || nrow(gseaResult@result) == 0) {
        writeReasonPng(
            file.path("GSEA_results/plots", paste0("GSEA_", category, ".png")),
            paste0("No significant GSEA results for ", category)
        )
        next
    }

    ## Save results to CSV
    write.csv(
        as.data.frame(gseaResult@result),
        file.path("GSEA_results/csv", paste0("GSEA_", category, ".csv")),
        row.names = FALSE
    )

    ## Prepare data for visualization
    gseaDf <- as_tibble(gseaResult@result) %>%
        mutate(phenotype = case_when(
            NES > 0 ~ "disease",
            NES < 0 ~ "healthy"
        ))

    gseaSigDf <- gseaDf %>%
        filter(abs(NES) > nesThreshold & p.adjust < padjGsea & pvalue < pvalThreshold) %>%
        arrange(p.adjust, decreasing = TRUE)

    ## Special handling for GO categories (C5)
    if (category == "C5") {

        goDf <- gseaSigDf %>%
            mutate(GO_type = case_when(
                grepl("^GOCC_", Description) ~ "CC",
                grepl("^GOMF_", Description) ~ "MF",
                grepl("^GOBP_", Description) ~ "BP",
                TRUE ~ "Other"
            )) %>%
            group_by(GO_type) %>%
            slice_max(order_by = p.adjust, n = ntopProcesses)

        for (goType in unique(goDf$GO_type)) {
            goTypeDf <- goDf %>%
                filter(GO_type == goType)

            ## Generate dot plot
            if (nrow(goTypeDf) > 0) {
                p <- dotplot(gseaResult, showCategory = goTypeDf$ID)
                ggsave(
                    file.path("GSEA_results/plots", paste0("GSEA_", category, "_", goType, "_dotplot.pdf")),
                    p + ggtitle(paste("GSEA Dot Plot -", category, ":", goType)),
                    width = 12,
                    height = 10
                )

                ## Generate GSEA enrichment plot
                if (length(goTypeDf$ID) >= 5) {
                    p <- gseaplot2(gseaResult, geneSetID = head(goTypeDf$ID, 5), pvalue_table = FALSE)
                    ggsave(
                        file.path("GSEA_results/plots", paste0("GSEA_", category, "_", goType, "_RESplot.pdf")),
                        p,
                        width = 12,
                        height = 10
                    )
                }

                ## Generate bubble plot
                p <- ggplot(goTypeDf, aes(x = phenotype, y = ID)) +
                    geom_point(aes(size = setSize, color = NES, alpha = -log10(p.adjust))) +
                    scale_color_gradient(low = "blue", high = "red") +
                    theme_bw() +
                    ggtitle(paste("GSEA Bubble Plot -", category, ":", goType))
                ggsave(
                    file.path("GSEA_results/plots", paste0("GSEA_", category, "_", goType, "_bubbleplot.pdf")),
                    p,
                    width = 12,
                    height = 10
                )
            }
        }

    } else {
        ## For non-GO categories

        ## Generate dot plot
        if (nrow(gseaSigDf) > 0) {
            p <- dotplot(gseaResult, showCategory = head(gseaSigDf$ID, ntopProcesses))
            ggsave(
                file.path("GSEA_results/plots", paste0("GSEA_", category, "_dotplot.pdf")),
                p + ggtitle(paste("GSEA Dot Plot -", category)),
                width = 12,
                height = 10
            )

            ## Generate GSEA enrichment plot
            if (nrow(gseaSigDf) >= 5) {
                p <- gseaplot2(gseaResult, geneSetID = head(gseaSigDf$ID, 5), pvalue_table = FALSE)
                ggsave(
                    file.path("GSEA_results/plots", paste0("GSEA_", category, "_RESplot.pdf")),
                    p,
                    width = 12,
                    height = 10
                )
            }

            ## Generate bubble plot
            p <- ggplot(head(gseaSigDf, ntopProcesses), aes(x = phenotype, y = ID)) +
                geom_point(aes(size = setSize, color = NES, alpha = -log10(p.adjust))) +
                scale_color_gradient(low = "blue", high = "red") +
                theme_bw() +
                ggtitle(paste("GSEA Bubble Plot -", category))
            ggsave(
                file.path("GSEA_results/plots", paste0("GSEA_", category, "_bubbleplot.pdf")),
                p,
                width = 12,
                height = 10
            )
        }
    }
}
