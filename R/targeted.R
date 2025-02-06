#' Run differential expression testing on a targeted subset of genes
#' @description This function runs differential expression analysis on a
#' subset of genes that are associated with specific Gene Ontology (GO) terms.
#'
#' @param go.terms A character vector containing the GO terms of interest
#' (e.g., c("GO:0000009", "GO:0042110")).
#' @param host.genes A matrix or data.frame containing host gene count data,
#' with rows as genes and columns as samples.
#' @param mtx A data.frame containing normalized microbial gene count data,
#' with rows as microbial genes and columns as samples.
#' @param mtx.features A character vector specifying the microbial gene
#' features to use.
#' @param padj A string specifying the p-value adjustment method.
#' Available options can be checked using `p.adjust.methods` (default: "fdr").
#' @param verbose Logical indicating whether to print progress messages
#' (default: FALSE).
#'
#' @return A data.frame with DESeq2 results for each microbial gene feature,
#' including columns for baseMean,
#' log2FoldChange, standard error (lfcSE),
#' statistical significance (stat), p-value, and adjusted p-value (padj).
#'
#' @export
#' @examples
#' set.seed(123)
#' GMM.matrix <- suppressMessages(get_GMM_matrix())
#' mtx.feature.names <- c("k02005_hly_d_family_secretion_protein",
#'                   "k02007_cobalt_nickel_transport_system_permease_protein",
#'                   "k00248_butyrate_gene_1",
#'                   "K00634_butyrate_gene_2",
#'                   "K00929_butyrate_gene_3")
#'
#' # simulate some host data
#' relevant.host.genes <- c("PIGV","ALG12")
#' nonrelevant.host.genes <- c("goofballgene", "sillygoosegene")
#' all.host.genes <- c(relevant.host.genes, nonrelevant.host.genes)
#' host.gene.counts <- data.frame(vapply(X = all.host.genes,
#'                                    FUN = function(x) {
#'                                          rnbinom(n = 1000,
#'                                          size = 400000, prob = runif(1))},
#'                                    FUN.VALUE = numeric(1000)
#'                                     )
#'                             )
#'
#' # simulate some mtx data
#' microbial.gene.counts <- data.frame(vapply(X=mtx.feature.names,
#'                                            FUN=function(x){
#'                                                rnorm(n=1000, mean=100)},
#'                                            FUN.VALUE=numeric(1000))
#'                                     )
#'
#' mtx.features.to.test <- features_from_gmm_df("butyrate production I",
#'                                        GMM.matrix,
#'                                        mtx.feature.names=mtx.feature.names)
#'
#' results <- GO_targeted_for_each_KO_within_GMM(go.terms=c("GO:0000009"),
#'                               host.genes=t(host.gene.counts),
#'                               mtx=microbial.gene.counts,
#'                               mtx.features=mtx.features.to.test)
#'
GO_targeted_for_each_KO_within_GMM <- function(go.terms, host.genes,
                                            mtx, mtx.features,
                                            padj="fdr",
                                            verbose=FALSE){
    # Get the host gene names to use
    full.targeted.genes <- c()
    for(go.term in go.terms){
        targeted.genes <- get_go_term_human_genes(go.term)
        # filter empty strings
        targeted.genes <- targeted.genes[nzchar(targeted.genes)]
        full.targeted.genes <- c(full.targeted.genes, targeted.genes)
    }
    full.targeted.genes <- unique(full.targeted.genes)

    # Initialize results object
    results <- data.frame(row.names=c("baseMean", "log2FoldChange", "lfcSE",
                                      "stat", "pvalue", "padj"))
    # This is kinda slow because it still relies on
    # pulling the GO terms each time
    for(feature in mtx.features){
        res <- go_targeted_diffex(targeted.genes=full.targeted.genes,
                                  host.genes=host.genes,
                                  microbial.gene=feature,
                                  microbial.genes=mtx,
                                  verbose=verbose)
        res$term <- feature
        results <- rbind(results, as.data.frame(res))

    }

    results$padj <- p.adjust(results$pvalue, method=padj)
    return(results)
}

#' Pull microbial features corresponding to a gut metabolic module (GMM)
#' @description Extracts all microbial transcriptomic (mtx) features
#' associated with a given gut metabolic module.
#'
#' @param GMM A string specifying the gut metabolic module of interest.
#' @param GMM.kos.df A data.frame with columns "Module" and "KEGG"
#' mapping KOs to their respective gut metabolic modules.
#' This can be generated using `get_GMM_matrix()`.
#' @param mtx.feature.names A character vector of all
#' metatranscriptomic feature names, including KEGG Orthology (KO) identifiers.
#'
#' @return A character vector of microbial features
#' corresponding to the specified GMM.
#' @export
#' @examples
#' GMM.matrix <- suppressMessages(get_GMM_matrix())
#' features <- features_from_gmm_df(GMM="butyrate production I",
#'                                  GMM.kos.df=GMM_matrix,
#'                                  mtx.feature.names=c("K00248", "K00634"))
#'
features_from_gmm_df <- function(GMM, GMM.kos.df, mtx.feature.names){
    module.kos <- GMM.kos.df[GMM.kos.df$Module==GMM, "KEGG"]

    mtx.features <- c()
    for(ko in module.kos){
        ko.feature <- mtx.feature.names[grepl(pattern=ko, x=mtx.feature.names,
                                              ignore.case=TRUE)]
        mtx.features <- c(mtx.features, ko.feature)
    }
    return(mtx.features)
}


#' Retrieve human genes associated with a GO term
#' @description Fetches all human genes that are annotated with a given
#' Gene Ontology (GO) term.
#'
#' @param go.term A character string representing a GO term,
#' such as "GO:0042110".
#'
#' @return A character vector of human gene symbols associated with the GO term.
#' @export
#'
#' @importFrom AnnotationDbi select
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#'
#' @examples
#' genes <- get_go_term_human_genes("GO:0000009")
#'
get_go_term_human_genes <- function(go.term){
    gene.data <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                       keys=c(go.term),
                                       columns=c('SYMBOL'),
                                       keytype="GOALL") # uses human
                                                        # ensembl annotations
    # gets gene symbol,
    # transcript_id and go_id for all genes annotated with the go term

    targeted.genes <- unique(gene.data$SYMBOL)
    return(targeted.genes)
}


#' Run targeted differential expression analysis using DESeq2
#' @description Performs differential expression analysis in DESeq2
#' on a subset of host genes,
#' based on their association with microbial gene expression.
#'
#' @param targeted.genes A character vector of host gene names
#' to include in the analysis.
#' @param host.genes A matrix or data.frame containing host gene count data,
#' with rows as genes and columns as samples.
#' @param microbial.gene A string specifying the microbial gene of interest.
#' @param microbial.genes A matrix or data.frame containing microbial
#' gene count data.
#' @param verbose Logical indicating whether to print progress messages
#' (default: TRUE).
#'
#' @return A data.frame with DESeq2 results,
#' including columns for baseMean, log2FoldChange,
#' standard error (lfcSE),
#' statistical significance (stat),
#' p-value, and adjusted p-value (padj).
#'
#' @export
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom DESeq2 DESeq
#' @importFrom DESeq2 results
#'
#' @examples
#' set.seed(123)
#' GMM.matrix <- suppressMessages(get_GMM_matrix())
#' mtx.feature.names <- c("k02005_hly_d_family_secretion_protein",
#'                   "k02007_cobalt_nickel_transport_system_permease_protein",
#'                   "k00248_butyrate_gene_1",
#'                   "K00634_butyrate_gene_2",
#'                   "K00929_butyrate_gene_3")
#'
#' # simulate some host data
#' relevant.host.genes <- c("PIGV","ALG12")
#' nonrelevant.host.genes <- c("goofballgene", "sillygoosegene")
#' all.host.genes <- c(relevant.host.genes, nonrelevant.host.genes)
#' host.gene.counts <- data.frame(vapply(X = all.host.genes,
#'                                    FUN = function(x) {
#'                                          rnbinom(n = 1000,
#'                                          size = 400000, prob = runif(1))},
#'                                    FUN.VALUE = numeric(1000)
#'                                     )
#'                             )
#'
#' # simulate some mtx data
#' microbial.gene <- "k00248_butyrate_gene_1"
#' microbial.gene.counts <- data.frame(vapply(X=mtx.feature.names,
#'                                            FUN=function(x){
#'                                                rnorm(n=1000, mean=100)},
#'                                            FUN.VALUE=numeric(1000))
#'                                     )
#'
#' results <- go_targeted_diffex(targeted.genes=relevant.host.genes,
#'                               host.genes=t(host.gene.counts),
#'                               microbial.gene=microbial.gene,
#'                               microbial.genes=microbial.gene.counts,
#'                               verbose=FALSE)
#'
go_targeted_diffex <- function(targeted.genes, host.genes,
                               microbial.gene, microbial.genes,
                               verbose=TRUE){
    targeted.data <- host.genes[targeted.genes,]
    if(verbose){
        s1 <- paste0(length(targeted.genes), " genes identified from GO term")
        s2 <- paste0(nrow(targeted.data), " genes in DESeq2 dataset")

        message(s1)
        message(s2)
    }

    ds <- DESeq2::DESeqDataSetFromMatrix(
        countData=targeted.data,
        colData=microbial.genes,
        design=as.formula(paste0("~",microbial.gene)))
    dds <- DESeq2::DESeq(ds)

    res <- DESeq2::results(dds, pAdjustMethod="fdr")
    return(res)
}

