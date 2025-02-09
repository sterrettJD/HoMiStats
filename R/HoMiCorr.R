#' Create a correlation network between host and microbial data
#' @description Runs regressions between all features
#' and returns the model outputs
#' @param mtx A dataframe of metatranscriptome gene counts,
#' where rows are samples, and columns are metatranscriptome genes.
#' Row names should be sample IDs.
#' @param host A dataframe of tpm-normalized host gene counts,
#' where rows are samples, and columns are host genes.
#' Row names should be sample IDs.
#' @param covariates A string formula with covariates to be included
#' in each regression.
#' All variables in this formula should be of a numeric type
#' (factors should be dummy coded, such that the reference value is 0).
#' @param metadata The metadata dataframe,
#' where rows are samples and columns are features.
#' Only needed if passing covariates.
#' @param sampleID A string denoting name of the column in metadata with
#' sample IDs,
#' which should be used to merge metadata with the feature table's rownames
#' @param reg.method A string denoting the method to use for regression.
#' Options include "zibr" (Zero-inflated beta regression with random effects),
#' "gamlss" (Zero-inflated beta regression implemented via GAMLSS),
#' "lm" (linear regression),
#' and "lmer" (linear mixed effects regression
#' implemented via lme4 and lmerTest).
#' @param padj A string denoting the p value adustment method.
#' Options can be checked using 'p.adjust.methods'
#' @param zero_prop_from_formula In ZIBR zero-inflated beta regression,
#' should the zeroes be modeled with the provided formula? Default is TRUE.
#' @param zibr_time_ind A string denoting the name of the time column
#' for ZIBR.
#' Defaults to NULL, which is implemented as a constant time value in ZIBR
#' to not fit a time effect.
#' This argument does nothing if reg.method is not "zibr".
#' @param ncores An integer denoting the number of cores to use if
#' running in parallel. Defaults to 1 (not parallelized).
#' @param show_progress A boolean denoting if a progress bar should be shown.
#' @return A dataframe with the model summaries as an adjacency list
#' @export
#' @importFrom Rfast comb_n
#' @importFrom broom.mixed tidy
#' @importFrom broom tidy
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @examples
#' # Simulate data
#' mtx <- data.frame(a1=c(0.1, 0.0, 0.0, 0.0),
#'                   b1=c(0.5, 0.5, 0.5, 0.4),
#'                   c1=c(0.4, 0.5, 0.0, 0.0),
#'                   d1=c(0.0, 0.0, 0.5, 0.6))
#' host <- data.frame(a2=c(0.5, 0.0, 0.0, 0.0),
#'                    b2=c(0.5, 0.5, 0.5, 0.4),
#'                    c2=c(0.4, 0.5, 0.0, 0.0),
#'                    d2=c(0.0, 0.0, 0.5, 0.6))
#' row.names(mtx) <- paste0("sample_", seq_len(4))
#' row.names(host) <- paste0("sample_", seq_len(4))
#' metadata <- data.frame(SampleID=paste0("sample_", seq_len(4)),
#'                        phenotype=c(0,0,1,1),
#'                        participant=c(0,1,0,1),
#'                        timepoint=c(0,0,1,1))
#'
#' # Run HoMiCorr using a zero inflated beta regression
#' # with random effects
#' homicorr.res.zibr <- run_HoMiCorr(mtx, host,
#'                                 reg.method="zibr",
#'                                 covariates="(1|participant)",
#'                                 metadata=metadata,
#'                                 sampleID="SampleID",
#'                                 zibr_time_ind="timepoint")
#'
#' # Run HoMiCorr using a (faster) zero inflated beta regression
#' # without random effects
#' homicorr.res.gamlss <- run_HoMiCorr(mtx, host,
#'                                 reg.method="gamlss",
#'                                 show_progress=FALSE)
#'
#'
run_HoMiCorr <- function(mtx, host, covariates=NULL, metadata=NULL,
                      sampleID=NULL, reg.method="zibr", padj="fdr",
                      zero_prop_from_formula=TRUE, zibr_time_ind=NULL,
                      ncores=1, show_progress=TRUE) {
    
    check_duplicated_colnames(colnames(mtx), colnames(host))
    if (reg.method %in% c("zibr", "gamlss")) {
        check_for_ones(mtx)
        check_for_ones(host)
    }
    
    data <- .prepare_data(mtx, host, covariates, metadata,
                        sampleID, zibr_time_ind)
    feature.combos <- .generate_feature_combinations(mtx, host)
    mod.summaries <- .initialize_model_summaries(feature.combos, reg.method)
    n.iterations <- length(feature.combos)

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl, cores=ncores)
    doSNOW::registerDoSNOW(cl)

    if (show_progress) {
        message("Running ", n.iterations, " iterations")
        pb <- txtProgressBar(max=n.iterations-1, style=3)
        doSNOWopts <- list(progress = function(n) setTxtProgressBar(pb, n))
    } else {
        pb <- NULL
        doSNOWopts <- list()
    }
    mod.summaries <- .run_regressions(feature.combos, data, reg.method, 
                                     covariates, zero_prop_from_formula, 
                                     zibr_time_ind, show_progress, 
                                     snowopts=doSNOWopts)
    .stop_parallel_processing(cl, show_progress, pb)
    
    return(.adjust_p_values(mod.summaries, reg.method,
                           zero_prop_from_formula, padj))
}


#' Prepare Data for HoMiCorr (internal)
#'
#' @description This function merges microbiome data (`mtx`),
#' host data (`host`), and optional metadata (`metadata`)
#' while selecting relevant covariates and sample identifiers.
#'
#' @param mtx A data frame containing microbiome data with row names as
#' sample identifiers and gene names as columns.
#' @param host A data frame containing host data with row names as
#' sample identifiers and gene names as columns.
#' @param covariates A character string specifying covariates to include
#' in the model formula (e.g., `"age + sex"`).
#' Random effects can also be passed (e.g., `"age + (1|participantID)"`).
#' If NULL, no covariates are used.
#' @param metadata A data frame containing metadata, including sample
#' identifiers and covariates.
#' Rows should correspond to samples, and covariates should be columns.
#' Can be NULL if no additional metadata is needed.
#' @param sampleID A character string specifying the column name in
#' `metadata` that contains sample identifiers.
#' @param zibr_time_ind A character string specifying a column in `metadata`
#' to be included in the merged data.
#'
#' @return A merged data frame containing microbiome data, host data,
#' and relevant metadata. The row names are set to sample identifiers.
#'
#' @examples
#' # Simulate data with samples as rows and genes as columns
#' mtx <- data.frame(a1 = c(0.1, 0.0, 0.0, 0.0),
#'                   b1 = c(0.5, 0.5, 0.5, 0.4))
#' rownames(mtx) <- paste0("sample_", seq_len(4))
#'
#' host <- data.frame(a2 = c(0.5, 0.0, 0.0, 0.0),
#'                    b2 = c(0.5, 0.5, 0.5, 0.4))
#' rownames(host) <- paste0("sample_", seq_len(4))
#'
#' metadata <- data.frame(SampleID = paste0("sample_", seq_len(4)),
#'                        age = c(25, 30, 35, 40),
#'                        sex = c("M", "F", "M", "F"),
#'                        participantID = c(0, 1, 0, 1),
#'                        timepoint = c(0, 0, 1, 1))
#' rownames(metadata) <- metadata$SampleID
#'
#' # Run .prepare_data function
#' result <- HoMiStats:::.prepare_data(mtx, host,
#'                                      "age + sex + (1|participantID)",
#'                                      metadata, "SampleID", "timepoint")
#' @keywords internal
#'
.prepare_data <- function(mtx, host,
                         covariates=NULL, metadata=NULL,
                         sampleID=NULL, zibr_time_ind=NULL) {
    
    # Check row names on mtx and host match each other
    if (length(setdiff(rownames(mtx), rownames(host))) > 0) {
        stop("mtx and host rownames do not match. ",
             "Please only make sure all samples in one dataset ",
             "exist in the other.")
    }
    formula <- if (!is.null(covariates)) paste0(" ~ ", covariates) else NULL
    metadata.vars <- c(all.vars(as.formula(formula)), sampleID)
    data <- merge(mtx, host, by="row.names")

    # Merge with metadata if applicable
    if ((!is.null(metadata)) && (!is.null(covariates))) {
        if ((sampleID %in% colnames(metadata)) == FALSE) {
            stop("sampleID column '", sampleID, "' not found in metadata.")
        }
        if (length(setdiff(rownames(mtx), metadata[, sampleID])) > 0) {
            stop("Row names of `mtx` do not match `metadata$", sampleID, "`.")
        }
        data <- merge(data, metadata[, c(metadata.vars, zibr_time_ind)],
                      by.x="Row.names", by.y=sampleID)
    }

    # Clean up output file
    row.names(data) <- data$Row.names
    data$Row.names <- NULL
    return(data)
}

#' Generate Feature Combinations for Regression Analysis (internal)
#'
#' @description This function generates all possible pairwise combinations
#' of microbial and host gene features.
#'
#' @param mtx A data frame containing microbial gene expression data,
#' where rows represent samples and columns represent genes.
#' @param host A data frame containing host gene expression data,
#' where rows represent samples and columns represent genes.
#'
#' @return A list of feature index pairs, representing all unique
#' combinations of microbial and host genes.
#'
#' @examples
#' # Simulated data
#' mtx <- data.frame(a1 = c(0.1, 0.0, 0.0, 0.0),
#'                   b1 = c(0.5, 0.5, 0.5, 0.4))
#' host <- data.frame(a2 = c(0.5, 0.0, 0.0, 0.0),
#'                    b2 = c(0.5, 0.5, 0.5, 0.4))
#'
#' # Generate feature combinations
#' feature_combos <- HoMiStats:::.generate_feature_combinations(mtx, host)
#'
#' @keywords internal⁠
#'
.generate_feature_combinations <- function(mtx, host) {
    all.featurenames <- c(colnames(mtx), colnames(host))
    return(Rfast::comb_n(seq_len(length(all.featurenames)),
                         k=2, simplify=FALSE))
}


#' Initialize Data Frame for Model Summaries (internal)
#'
#' @description Creates an empty data frame to store regression model summaries.
#' The columns included in the data frame depend on the selected
#' regression method.
#'
#' @param feature.combos A list of feature combinations to be analyzed.
#' @param reg.method A string specifying the regression method to be used.
#' Supported options: "gamlss", "zibr", "lm", "lmer".
#'
#' @return A data frame initialized with appropriate columns for the chosen
#' regression method, with the number of rows equal to the number of feature
#' combinations.
#'
#' @examples
#' # Simulated feature combinations
#' feature_combos <- list(c("a1", "b2"), c("b1", "a2"))
#'
#' # Initialize model summaries for linear regression
#' summaries <- HoMiStats:::.initialize_model_summaries(feature_combos, "lm")
#'
#' @keywords internal⁠
#'
.initialize_model_summaries <- function(feature.combos, reg.method) {
    n.iterations <- length(feature.combos)
    summary.cols <- switch(reg.method,
        "gamlss" = c("parameter", "term", "estimate",
                     "std.error", "statistic", "p.value", "feature"),
        "zibr" = c("parameter", "term", "estimate",
                   "p.value", "joint.p", "feature"),
        "lm" = c("term", "estimate", "std.error",
                 "statistic", "p.value", "feature"),
        "lmer" = c("effect", "group", "term", "estimate",
                   "std.error", "statistic", "df", "p.value", "feature"))
    
    return(data.frame(matrix(nrow=n.iterations, 
                            ncol=length(summary.cols), 
                            dimnames=list(NULL, summary.cols))))
}
#' Stop Parallel Processing and Clean Up (internal)
#'
#' @description This function stops the parallel processing cluster
#' and closes the progress bar if `show_progress` is enabled.
#'
#' @param cl A parallel cluster object created by `parallel::makeCluster`.
#' @param show_progress A boolean indicating whether a progress bar was shown.
#' If TRUE, the function attempts to close the progress bar.
#'
#' @return This function does not return a value. It stops the cluster
#' and cleans up resources.
#'
#' @examples
#' # Example usage (not run):
#' # cl <- parallel::makeCluster(2)
#' # HoMiStats:::.stop_parallel_processing(cl, show_progress = FALSE)
#'
#' @keywords internal
#'
.stop_parallel_processing <- function(cl, show_progress, pb) {
    if (show_progress) close(pb)
    parallel::stopCluster(cl)
}


#' Adjust P-values for Multiple Comparisons (internal)
#'
#' @description This function applies multiple testing correction
#' to the p-values from regression model summaries. The method used for
#' adjustment is specified by the user.
#'
#' @param mod.summaries A data frame containing model summary results,
#' including p-values.
#' @param reg.method A string indicating the regression method used.
#' Options include "zibr", "gamlss", "lm", and "lmer".
#' @param zero_prop_from_formula A boolean indicating whether zeroes
#' were modeled using the provided formula (relevant for ZIBR regression).
#' @param padj A string specifying the multiple testing correction method.
#' Supported methods can be checked with `p.adjust.methods`.
#'
#' @return A modified version of `mod.summaries` with an additional column
#' `q`, containing the adjusted p-values.
#'
#' @examples
#' # Example model summary
#' mod.summaries <- data.frame(term = c("(Intercept)", "GeneX"),
#'                             parameter = c("beta", "beta"),
#'                             p.value = c(0.01, 0.002),
#'                             joint.p = c(0.02, 0.005))
#'
#' # Adjust p-values using FDR correction
#' mod.summaries <- HoMiStats:::.adjust_p_values(mod.summaries, "zibr",
#'                                  zero_prop_from_formula = TRUE,
#'                                  padj = "fdr")
#'
#' @keywords internal
#'
.adjust_p_values <- function(mod.summaries, reg.method,
                            zero_prop_from_formula, padj) {
    mod.summaries <- as.data.frame(mod.summaries)
    # adjust p value only for non-intercept terms
    # and if we have a joint p, only adjust it for the beta coefficient
    # (it's copied for the logistic)
    if((reg.method == "zibr") & zero_prop_from_formula){
        mod.summaries[mod.summaries$term!="(Intercept)" &
                      mod.summaries$parameter == "beta",
                      "q"] <- p.adjust(
                          mod.summaries[mod.summaries$term!="(Intercept)" &
                          mod.summaries$parameter == "beta",
                                        "joint.p"],
                          method=padj)

    } else if((reg.method == "gamlss") |
              (reg.method == "zibr") & (zero_prop_from_formula == FALSE)){
        mod.summaries[mod.summaries$term!="(Intercept)",
                      "q"] <- p.adjust(
                          mod.summaries[mod.summaries$term!="(Intercept)",
                                        "p.value"],
                          method=padj)

    } else if(reg.method == "lmer"){
        mod.summaries[mod.summaries$term!="(Intercept)" &
                      mod.summaries$effect == "fixed",
                      "q"] <- p.adjust(
                          mod.summaries[mod.summaries$term!="(Intercept)" &
                          mod.summaries$effect == "fixed",
                                        "p.value"],
                          method=padj)

    } else if(reg.method == "lm"){
        mod.summaries[mod.summaries$term!="(Intercept)",
                      "q"] <- p.adjust(
                          mod.summaries[mod.summaries$term!="(Intercept)",
                                        "p.value"],
                          method=padj)
    }
    return(mod.summaries)
}

#' Run Pairwise Regressions Across Features (internal)
#'
#' @description This function iterates over all specified feature pairs
#' and performs regression analyses in parallel. The regression method
#' is determined by `reg.method`, and the function supports multiple
#' modeling approaches, including ZIBR, GAMLSS, linear models,
#' and linear mixed-effects models.
#'
#' @param feature.combos A list of feature index pairs to be analyzed.
#' Each pair specifies two columns in `data` to be used as dependent and
#' independent variables.
#' @param data A data frame containing the feature columns and covariates
#' for modeling.
#' @param reg.method A string specifying the regression method to use.
#' Options include `"zibr"`, `"gamlss"`, `"lm"`, and `"lmer"`.
#' @param covariates A string specifying the covariates to include
#' in the model formula.
#' @param zero_prop_from_formula A boolean indicating whether
#' zero proportions should be modeled
#' using the specified formula (relevant for ZIBR regression).
#' @param zibr_time_ind A string specifying the time variable for ZIBR models.
#' @param show_progress A boolean indicating whether to
#' display progress updates.
#' @param snowopts Options passed to `foreach::foreach()`
#' for parallel execution using the `doSNOW` backend.
#'
#' @return A data frame summarizing the regression results
#' for each feature pair.
#'
#' @examples
#' # Example usage (not run):
#' # feature.combos <- list(c(1, 2), c(3, 4))
#' # results <- HoMiStats:::.run_regressions(feature.combos, data,
#' #                             "lm", "age + sex",
#' #                             NULL, NULL, TRUE)
#'
#' @keywords internal
#'
.run_regressions <- function(feature.combos, data, reg.method,
                            covariates, zero_prop_from_formula,
                            zibr_time_ind, show_progress, snowopts=doSNOWopts) {
    mod.summaries <- foreach::foreach(cols=feature.combos, .combine=rbind,
        .options.snow=snowopts,
        .export = c(".run_single_model")) %dopar% {
        col1 <- colnames(data)[cols[1]]
        col2 <- colnames(data)[cols[2]]
        mod.sum <- .run_single_model(col1, col2, data, reg.method, covariates,
                                    zero_prop_from_formula, zibr_time_ind)
        if(nrow(mod.sum) > 0) {
            mod.sum$feature <- col1
        }
        return(mod.sum)
    }
    
    return(mod.summaries)
}



#' Run a Single Regression Model
#'
#' @description This function performs a single regression analysis
#' for a given pair of features using the specified regression method.
#'
#' @param col1 A string specifying the dependent variable (feature) in `data`.
#' @param col2 A string specifying the independent variable (feature) in `data`.
#' @param data A data frame containing the feature columns and covariates.
#' @param reg.method A string indicating the type of regression to run.
#' Options include `"gamlss"`, `"zibr"`, `"lm"`, and `"lmer"`.
#' @param covariates A string specifying the covariates to
#' include in the model formula.
#' @param zero_prop_from_formula A boolean indicating whether
#' zero proportions should be modeled using the provided formula
#' (relevant for ZIBR regression).
#' @param zibr_time_ind A string specifying the time variable for ZIBR models.
#'
#' @return A data frame containing the regression summary for `col2`,
#' including coefficient estimates, standard errors, and p-values.
#' If the model fails to converge (for mixed-effects models), an error message
#' is returned instead.
#'
#' @examples
#' # Example usage (not run):
#' # result <- HoMiStats:::.run_single_model("geneA", "geneB", data, "lm",
#' #                            "age + sex", FALSE, NULL)
#'
#' @keywords internal
#'
.run_single_model <- function(col1, col2, data, reg.method, covariates,
                             zero_prop_from_formula, zibr_time_ind) {
    if (reg.method == "gamlss") {
        mod <- run_single_beta_reg_gamlss(paste0(col1, " ~ ", 
                                                covariates, " + ", col2), data)
        mod.sum <- broom.mixed::tidy(mod)
        return(mod.sum[mod.sum$term == col2 & mod.sum$parameter == "mu", ])
    }
    if (reg.method == "zibr") {
        if(!is.null(covariates)) {
            form <- as.formula(paste0("~ ", covariates))
            random.effects.vars <- get_random_fx(form)
            fixed.vars <- setdiff(all.vars(form), random.effects.vars)
        } else {
            random.effects.vars <- NULL
            fixed.vars <- NULL
        }
        log.cov <- if (zero_prop_from_formula) c(fixed.vars, col2) else NULL
        mod <- run_single_beta_reg_zibr(logistic_cov=log.cov,
                                        beta_cov=c(fixed.vars, col2),
                                        Y=col1,
                                        subject_ind=random.effects.vars,
                                        time_ind=zibr_time_ind,
                                        data=data)

        mod.sum <- tidy_zibr_results(mod)
        mod.sum$term <- map_zibr_termnames(mod.sum$term, c(fixed.vars, col2))
        return(mod.sum[mod.sum$term == col2 & mod.sum$parameter == "beta", ])
    }
    if (reg.method == "lm") {
        mod <- run_single_lm(paste0(col1, " ~ ", covariates, " + ", col2), data)
        mod.sum <- broom::tidy(mod)
        return(mod.sum[mod.sum$term == col2, ])
    }
    if (reg.method == "lmer") {
        # This try catch block is used for convergence (or similar) errors
        # But maybe should be more specific
        mod.sum <- tryCatch({
            mod <- run_single_lmer(paste0(col1, " ~ ", covariates, " + ", col2),
                                   data)
            broom.mixed::tidy(mod)
          }, error=function(e) {
            data.frame("effect"="error", "group"=NA,
                       "term"=col2, "estimate"=NA,
                       "std.error"=NA, "statistic"=NA,
                       "df"=NA, "p.value"=NA)}) # END TRYCATCH
        # grab only the col2 beta row
        return(mod.sum[mod.sum$term==col2 & mod.sum$effect=="fixed",])
    }
}



#' Check for duplicated values between 2 sets of column names
#' @description Ensures there are no intersecting values between two vectors
#' @param colnames1 vector of strings with column names from the first dataset
#' @param colnames2 vector of strings with column names from the second dataset
#' @return Nothing
#' @export
#' @examples
#' # Simulate data
#' mtx <- data.frame(a1=c(0.1, 0.0, 0.0, 0.0),
#'                   b1=c(0.5, 0.5, 0.5, 0.4),
#'                   c1=c(0.4, 0.5, 0.0, 0.0),
#'                   d1=c(0.0, 0.0, 0.5, 0.6))
#' host <- data.frame(a2=c(0.5, 0.0, 0.0, 0.0),
#'                    b2=c(0.5, 0.5, 0.5, 0.4),
#'                    c2=c(0.4, 0.5, 0.0, 0.0),
#'                    d2=c(0.0, 0.0, 0.5, 0.6))
#'
#' # This should not raise an error
#' check_duplicated_colnames(colnames(mtx), colnames(host))
#' # This should raise an error
#' try(check_duplicated_colnames(colnames(mtx), c("a1")))
#'
check_duplicated_colnames <- function(colnames1, colnames2){
    colnames.intersection <- intersect(colnames1, colnames2)
    if(length(colnames.intersection) > 0){
        colnames.intersection.chr <- paste(colnames.intersection, collapse=", ")
        stop("Columns found in both datasets: ", colnames.intersection.chr)
    }
}
