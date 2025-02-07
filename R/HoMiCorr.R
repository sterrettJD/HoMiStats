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
#' @param zibr_zibr_time_ind A string denoting the name of the time column
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
    
    data <- prepare_data(mtx, host, covariates, metadata, sampleID, zibr_time_ind)
    feature.combos <- generate_feature_combinations(mtx, host)
    mod.summaries <- initialize_model_summaries(feature.combos, reg.method)
    n.iterations <- length(feature.combos)

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl, cores=ncores)
    doSNOW::registerDoSNOW(cl)

    if (show_progress) {
        message("Running ", n.iterations, " iterations")
        pb <- txtProgressBar(max=n.iterations-1, style=3)
        doSNOWopts <- list(progress = function(n) setTxtProgressBar(pb, n))
    } else {
        doSNOWopts <- list()
    }
    mod.summaries <- run_regressions(feature.combos, data, reg.method, 
                                     covariates, zero_prop_from_formula, 
                                     zibr_time_ind, show_progress, 
                                     snowopts=doSNOWopts)
    stop_parallel_processing(cl, show_progress)
    
    return(adjust_p_values(mod.summaries, reg.method,
                           zero_prop_from_formula, padj))
}

prepare_data <- function(mtx, host, covariates, metadata, sampleID, zibr_time_ind) {
    formula <- if (!is.null(covariates)) paste0(" ~ ", covariates) else NULL
    metadata.vars <- c(all.vars(as.formula(formula)), sampleID)
    data <- merge(mtx, host, by="row.names")
    if (!is.null(metadata)) {
        data <- merge(data, metadata[, c(metadata.vars, zibr_time_ind)], 
                      by.x="Row.names", by.y=sampleID)
    }
    row.names(data) <- data$Row.names
    data$Row.names <- NULL
    return(data)
}

generate_feature_combinations <- function(mtx, host) {
    all.featurenames <- c(colnames(mtx), colnames(host))
    return(Rfast::comb_n(seq_len(length(all.featurenames)),
                         k=2, simplify=FALSE))
}

initialize_model_summaries <- function(feature.combos, reg.method) {
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

stop_parallel_processing <- function(cl, show_progress) {
    if (show_progress) close(pb)
    parallel::stopCluster(cl)
}

adjust_p_values <- function(mod.summaries, reg.method, zero_prop_from_formula, padj) {
    term_filter <- mod.summaries$term != "(Intercept)"
    if (reg.method == "zibr" & zero_prop_from_formula) {
        term_filter <- term_filter & mod.summaries$parameter == "beta"
        mod.summaries$q <- p.adjust(mod.summaries$joint.p[term_filter], method=padj)
    } else {
        mod.summaries$q <- p.adjust(mod.summaries$p.value[term_filter], method=padj)
    }

    return(mod.summaries)
}

run_regressions <- function(feature.combos, data, reg.method, covariates, zero_prop_from_formula, zibr_time_ind, show_progress, snowopts=doSNOWopts) {
    mod.summaries <- foreach::foreach(cols=feature.combos, .combine=rbind, 
        .options.snow=snowopts,
        .export = c("run_single_model")) %dopar% {
        col1 <- colnames(data)[cols[1]]
        col2 <- colnames(data)[cols[2]]
        mod.sum <- run_single_model(col1, col2, data, reg.method, covariates, zero_prop_from_formula, zibr_time_ind)
        
        if(nrow(mod.sum) > 0) {mod.sum$feature <- col1}
        
        return(mod.sum)
    }
    
    return(mod.summaries)
}

run_single_model <- function(col1, col2, data, reg.method, covariates, zero_prop_from_formula, zibr_time_ind) {
    if (reg.method == "gamlss") {
        mod <- run_single_beta_reg_gamlss(paste0(col1, " ~ ", covariates, " + ", col2), data)
        return(broom.mixed::tidy(mod)[broom.mixed::tidy(mod)$term == col2 & broom.mixed::tidy(mod)$parameter == "mu", ])
    }
    if (reg.method == "zibr") {
        if(!is.null(covariates)){
            # extract vars for zibr
            form <- as.formula(paste0("~ ", covariates))
            vars <- all.vars(form)
            # Extracts random effects from formula
            random.effects.vars <- get_random_fx(form)
            fixed.vars <- setdiff(vars, random.effects.vars)
        } else {
            # Extracts random effects from formula
            random.effects.vars <- NULL
            fixed.vars <- NULL
        }

        mod <- run_single_beta_reg_zibr(logistic_cov=if (zero_prop_from_formula) c(fixed.vars, col2) else NULL,
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
        mod.sum <- tryCatch({
            mod <- run_single_lmer(paste0(col1, " ~ ", covariates, " + ", col2),
                                   data)
            broom.mixed::tidy(mod)
          }, error=function(e) {
            data.frame("effect"="error", "group"=NA,
                                    "term"=col2, "estimate"=NA,
                                    "std.error"=NA, "statistic"=NA,
                                    "df"=NA, "p.value"=NA)})
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
#' check_duplicated_colnames(colnames(mtx), c("a1"))
#'
check_duplicated_colnames <- function(colnames1, colnames2){
    colnames.intersection <- intersect(colnames1, colnames2)
    if(length(colnames.intersection) > 0){
        colnames.intersection.chr <- paste(colnames.intersection, collapse=", ")
        stop("Columns found in both datasets: ", colnames.intersection.chr)
    }
}


