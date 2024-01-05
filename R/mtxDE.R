library(gamlss)

# runs a zero-inflated beta distribution GAM using gamlss
# returns the model

#' Run a single zero-inflatd beta regression using gamlss
#' @description runs a zero-inflated beta distribution GAM using gamlss
#' @param formula the full formula to be used in the regression
#' @param data A dataframe to be used in the regression
#' @return the resulting model
#' @export
#' @importFrom gamlss.dist BEZI
#'
run_single_beta_regression <- function(formula, data){
    mod <- gamlss::gamlss(as.formula(formula),
                        data=data,
                        family=gamlss.dist::BEZI)
    return(mod)
}

#' Run differential expression analysis for metatranscriptomics data
#' @description Regresses each feature using the formula provided and returns the statistical summary for each feature
#' @param formula the right hand side of the regression formula to be used for each feature
#' @param feature.table A dataframe, where rows are samples, and columns are genes/features. Row names should be sample IDs.
#' @param metadata The metadata dataframe, where rows are samples and columns are features
#' @param sampleID The name of the column in metadata with sample IDs, which should be used to merge metadata with the feature table's rownames
#' @param padj p value adustment method. Options can be checked using 'p.adjust.methods'
#' @return a dataframe with differential expression results
#' @export
#' @importFrom broom.mixed tidy
#'
run_mtxDE <- function(formula, feature.table, metadata, sampleID, padj){
    # merge the feature table and metadata based on the rownames
    formula <- paste0(" ~ ", formula)
    metadata.vars <- c(all.vars(as.formula(formula)), sampleID)
    data <- merge(feature.table, metadata[,metadata.vars], by="row.names")

    # initialize the model summaries df
    mod.summaries <- data.frame(matrix(nrow=0, ncol=7))
    colnames(mod.summaries) <- c("parameter", "term",
                                 "estimate", "std.error",
                                 "statistic", "p.value",
                                 "feature")
    # Loop through each column and run the beta regression
    for(col in colnames(feature.table)){
        mod <- run_single_beta_regression(paste0(col, formula),
                                          data=data)
        mod.sum <- broom.mixed::tidy(mod)
        mod.sum$feature <- col
        mod.summaries <- rbind(mod.summaries, mod.sum)
    }
    mod.summaries$q <- p.adjust(mod.summaries$p.value, method=padj)
    return(mod.summaries)
}
