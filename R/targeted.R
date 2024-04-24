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
GO_targeted_for_each_KO_within_GMM <- function(go.term, host.genes,
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

#' Pull all human genes corresponding to a GO term of interest
#' @description Looks for all genes corresponding to the GO term of interest
#' @param go.term A string version of a GO term, such as "GO:0042110"
#' @return a vector of gene names corresponding to the GO term of interest
#' @export
#' @importFrom AnnotationDbi select
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#'
get_go_term_human_genes <- function(go.term){
    gene.data <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                       keys=c(go.term),
                                       columns=c('SYMBOL'),
                                       keytype="GOALL") #uses human ensembl annotations
    #gets gene symbol, transcript_id and go_id for all genes annotated with the go term

    targeted.genes <- unique(gene.data$SYMBOL)
    return(targeted.genes)
}


#' Run targeted differential expression analysis in DESeq2 for a subset of genes
#' @description Uses the vector of genes provided to subset the data then run DESeq2
#' @param targeted.genes A vector of host gene names
#' @param host.genes A matrix or data.frame with the host gene count data
#' @param microbial.gene The name of the microbial gene to be used
#' @param microbial.genes A matrix or data.frame with the microbial gene count data
#' @return A dataframe of DESeq2 results
#' @export
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom DESeq2 DESeq2
#' @importFrom DESeq2 results
#'
go_targeted_diffex <- function(targeted.genes, host.genes,
                               microbial.gene, microbial.genes,
                               verbose=T){
    targeted.genes <- get_go_term_human_genes(go.term)

    # filter empty strings
    targeted.genes <- targeted.genes[nzchar(targeted.genes)]


    targeted.data <- host.genes[targeted.genes,]
    if(verbose){
        s1 <- paste0(length(targeted.genes), " genes identified from GO term")
        s2 <- paste0(nrow(targeted.data), " genes in DESeq2 dataset")
        print(s1)
        print(s2)
    }

    ds <- DESeq2::DESeqDataSetFromMatrix(
        countData=targeted.data,
        colData=microbial.genes,
        design=as.formula(paste0("~",microbial.gene)))
    dds <- DESeq2::DESeq(ds)

    res <- DESeq2::results(dds, pAdjustMethod="fdr")
    return(res)
}

