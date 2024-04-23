#' Run differential expression testing on a targeted subset of genes
#' @description
#' @param go.term the GO term of interest
#' @param host.genes A matrix or data.frame with the host gene count data
#' @param mtx A data.frame with the normalized microbial gene count data
#' @param mtx.features A vector of features to use from the microbial gene counts
#' @param padj A string denoting the p value adustment method. Options can be checked using 'p.adjust.methods'
#' @return a data.frame with DESeq results for each mtx feature
#' @export
#'
targeted_for_each_KO <- function(go.term, host.genes,
                                 mtx, mtx.features,
                                 padj="fdr"){
    results <- data.frame(row.names=c("baseMean", "log2FoldChange", "lfcSE",
                                      "stat", "pvalue", "padj"))
    # This is kinda slow because it still relies on pulling the GO terms each time
    for(feature in mtx.features){
        res <- go_targeted_diffex(go.term, host.genes=gns.only.subset,
                                  microbial.genes=mtx,
                                  microbial.gene=feature)
        res$term <- feature
        results <- rbind(results, as.data.frame(res))

    }
    results$padj <- p.adjust(results$pvalue, method=padj)
    return(results)
}

#' Pull mtx features for a gut metabolic module of interest
#' @description Looks for all mtx features corresponding to the gut metabolic module of interest
#' @param GMM The name of the gut metabolic module of interest
#' @param GMM.kos.df A data.frame with columns "Module" and "KEGG" with the KOs corresponding to each gut metabolic module. Can be generated using get_GMM_matrix()
#' @param mtx A vector with all of the metatranscriptome feature names, with KOs in them
#' @return a vector of feature names corresponding to the GMM of interest
#' @export
#'
features_from_gmm_df <- function(GMM, GMM.kos.df, mtx.feature.names){
    module.kos <- GMMs.kos.df[GMMs.kos.df$Module==GMM, "KEGG"]

    mtx.features <- c()
    for(ko in module.kos){
        ko.feature <- mtx.feature.names[grepl(pattern=ko, x=mtx.feature.names,
                                              ignore.case=TRUE)]
        mtx.features <- c(mtx.features, ko.feature)
    }
    return(mtx.features)
}
