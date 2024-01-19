library(gamlss)

#' Run a single zero-inflatd beta regression using gamlss
#' @description runs a zero-inflated beta distribution GAM using gamlss
#' @param formula the full formula to be used in the regression
#' @param data A dataframe to be used in the regression
#' @return the resulting model
#' @export
#' @importFrom gamlss.dist BEZI
#'
run_single_beta_regression_gamlss <- function(formula, data,
                                              controller=gamlss.control(trace=FALSE)){
    mod <- gamlss::gamlss(as.formula(formula),
                        data=data,
                        family=gamlss.dist::BEZI,
                        control=controller)
    return(mod)
}

#' Check for 1 in feature.table
#' @description The zero-inflated beta regression can't handle values of 1. This checks for them and raises an error if they exist.
#' @param feature.table A dataframe, where rows are samples, and columns are genes/features. Row names should be sample IDs.
#' @return Nothing
#' @export
#'
check_for_ones <- function(feature.table){
    rows.with.ones <- which(feature.table == 1, arr.ind = T)[,1]
    unique.rows.with.ones <- unique(rows.with.ones)
    unique.rownames.with.ones <- rownames(feature.table)[unique.rows.with.ones]
    if(length(unique.rownames.with.ones) > 0){
        stop(paste("The following rows contains a value of one which the zero-inflated beta regression cannot handle:",
                   unique.rownames.with.ones))
    }
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
run_mtxDE <- function(formula, feature.table, metadata, sampleID, padj="fdr"){
    # Check for values of one, which the beta regression can't handle
    check_for_ones(feature.table)

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

    n_iter <- ncol(feature.table)
    i <- 0
    # Initializes progress bar
    pb <- txtProgressBar(min = 0,
                         max = n_iter,
                         style = 3,
                         width = 50,   # Progress bar width. Defaults to getOption("width")
                         char = "=")

    # Loop through each column and run the beta regression
    for(col in colnames(feature.table)){
        mod <- run_single_beta_regression_gamlss(paste0(col, formula),
                                          data=data)
        mod.sum <- broom.mixed::tidy(mod)
        mod.sum$feature <- col
        mod.summaries <- rbind(mod.summaries, mod.sum)

        # progress bar things
        setTxtProgressBar(pb, i)
        i <- i + 1
    }

    close(pb) #close the progress bar

    mod.summaries <- as.data.frame(mod.summaries)
    # adjust p value only for non-intercept terms
    mod.summaries[mod.summaries$term!="(Intercept)",
                  "q"] <- p.adjust(mod.summaries[mod.summaries$term!="(Intercept)",
                                                "p.value"],
                                   method=padj)

    return(mod.summaries)
}
