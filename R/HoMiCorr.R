#' Create a correlation network between host and microbial data
#' @description Runs regressions between all features and returns the model outputs
#' @param mtx A dataframe of metatranscriptome gene counts, where rows are samples, and columns are metatranscriptome genes. Row names should be sample IDs.
#' @param host A dataframe of tpm-normalized host gene counts, where rows are samples, and columns are host genes. Row names should be sample IDs.
#' @param covariates A string formula with covariates to be included in each regression. All variables in this formula should be of a numeric type (factors should be dummy coded, such that the reference value is 0).
#' @param metadata The metadata dataframe, where rows are samples and columns are features. Only needed if
#' @param sampleID A string denoting name of the column in metadata with sample IDs, which should be used to merge metadata with the feature table's rownames
#' @param reg.method A string denoting the method to use for regression. Options include "zibr" (Zero-inflated beta regression with random effects), "gamlss" (Zero-inflated beta regression implemented via GAMLSS), "lm" (linear regression), and "lmer" (linear mixed effects regression implemented via lme4 and lmerTest).
#' @param padj A string denoting the p value adustment method. Options can be checked using 'p.adjust.methods'
#' @param zero_prop_from_formula In ZIBR zero-inflated beta regression, should the zeroes be modeled with the provided formula? Default is TRUE.
#' @param zibr_zibr_time_ind A string denoting the name of the time column for ZIBR. Defaults to NULL, which is implemented as a constant time value in ZIBR to not fit a time effect. This argument does nothing if reg.method is not "zibr".
#' @param ncores An integer denoting the number of cores to use if running in parallel. Defaults to 1 (not parallelized).
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
#' row.names(mtx) <- paste0("sample_", 1:4)
#' row.names(host) <- paste0("sample_", 1:4)
#' metadata <- data.frame(SampleID=paste0("sample_", 1:4),
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
run_HoMiCorr <- function(mtx, host,
                      covariates=NULL, metadata=NULL,
                      sampleID=NULL,
                      reg.method="zibr",
                      padj="fdr",
                      zero_prop_from_formula=TRUE,
                      zibr_time_ind=NULL,
                      ncores=1,
                      show_progress=TRUE){

    # Make sure there aren't duplicated column names between the two datasets
    check_duplicated_colnames(colnames(mtx), colnames(host))
    # Check for values of one, which the beta regression can't handle
    if(reg.method %in% c("zibr", "gamlss")){
        check_for_ones(mtx)
        check_for_ones(host)
    }

    # merge the feature table and metadata based on the rownames
    if(!is.null(covariates)){
        formula <- paste0(" ~ ", covariates)
    } else {
        formula <- NULL
    }

    metadata.vars <- c(all.vars(as.formula(formula)), sampleID)

    data <- merge(mtx, host, by="row.names")
    if(!is.null(metadata)){
        data <- merge(data,
                      metadata[,c(metadata.vars, zibr_time_ind)],
                      by.x="Row.names", by.y=sampleID)
    }

    # calculate regressions to run
    all.featurenames <- c(colnames(mtx), colnames(host))
    n.feat <- length(all.featurenames)
    feature.combos <- Rfast::comb_n(1:n.feat,
                                    k=2,
                                    simplify=FALSE)
    n.iterations <- length(feature.combos)

    # initialize the model summaries df
    if(reg.method == "gamlss"){
        mod.summaries <- data.frame(matrix(nrow=n.iterations, ncol=7))
        colnames(mod.summaries) <- c("parameter", "term",
                                     "estimate", "std.error",
                                     "statistic", "p.value",
                                     "feature")

    } else if(reg.method == "zibr"){
        mod.summaries <- data.frame(matrix(nrow=n.iterations, ncol=6))
        colnames(mod.summaries) <- c("parameter", "term",
                                     "estimate",
                                     "p.value", "joint.p",
                                     "feature")
        if(!is.null(covariates)){
            # extract vars for zibr
            vars <- all.vars(as.formula(formula))
            # Extracts random effects from formula
            random.effects.vars <- get_random_fx(as.formula(formula))
            fixed.vars <- setdiff(vars, random.effects.vars)
        } else {
            # extract vars for zibr
            vars <- NULL
            # Extracts random effects from formula
            random.effects.vars <- NULL
            fixed.vars <- NULL
        }


    } else if(reg.method=="lm"){
        mod.summaries <- data.frame(matrix(nrow=n.iterations, ncol=6))
        colnames(mod.summaries) <- c("term", "estimate",
                                     "std.error", "statistic",
                                     "p.value", "feature")

    } else if(reg.method=="lmer"){
        mod.summaries <- data.frame(matrix(nrow=n.iterations, ncol=9))
        colnames(mod.summaries) <- c("effect", "group",
                                     "term", "estimate",
                                     "std.error", "statistic",
                                     "df", "p.value",
                                     "feature")

    }

    # Setup cluster for parallel running
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl, cores=ncores)
    doSNOW::registerDoSNOW(cl)

    if(show_progress){
        print(paste("Running", n.iterations, "interations"))
        # Initializes progress bar

        pb <- txtProgressBar(max=n.iterations-1, style=3)
        progress <- function(n) setTxtProgressBar(pb, n)
        doSNOWopts <- list(progress = progress)
    } else {
        doSNOWopts <- list()
    }

    # Loop through each combination of columns and run the regression
    mod.summaries <- foreach::foreach(cols=feature.combos,
                                      .combine=rbind,
                                      .options.snow=doSNOWopts) %dopar% {
        col1 <- all.featurenames[cols[1]]
        col2 <- all.featurenames[cols[2]]

      if(reg.method == "gamlss"){
          mod <- run_single_beta_reg_gamlss(paste0(col1, "~", covariates, " + ", col2),
                                            data=data)
          mod.sum <- broom.mixed::tidy(mod)
          # grab only the col2 beta row
          mod.sum <- mod.sum[mod.sum$term==col2 & mod.sum$parameter=="mu",]

      } else if((reg.method == "zibr") & (zero_prop_from_formula==TRUE)){
          mod <- run_single_beta_reg_zibr(logistic_cov=c(fixed.vars, col2),
                                          beta_cov=c(fixed.vars, col2),
                                          Y=col1,
                                          subject_ind=random.effects.vars,
                                          time_ind=zibr_time_ind,
                                          data=data)

          mod.sum <- tidy_zibr_results(mod)
          mod.sum$term <- map_zibr_termnames(mod.sum$term, c(fixed.vars, col2))
          # grab only the col2 beta row
          mod.sum <- mod.sum[mod.sum$term==col2 & mod.sum$parameter=="beta",]

      } else if((reg.method == "zibr") & (zero_prop_from_formula==FALSE)){

          mod <- run_single_beta_reg_zibr(logistic_cov=NULL,
                                          beta_cov=c(fixed.vars, col2),
                                          Y=col1,
                                          subject_ind=random.effects.vars,
                                          time_ind=zibr_time_ind,
                                          data=data)

          mod.sum <- tidy_zibr_results(mod)
          mod.sum$term <- map_zibr_termnames(mod.sum$term, c(fixed.vars, col2))
          # grab only the col2 beta row
          mod.sum <- mod.sum[mod.sum$term==col2 & mod.sum$parameter=="beta",]

      } else if(reg.method == "lm"){
          mod <- run_single_lm(paste0(col1, "~", covariates, " + ", col2), data)
          mod.sum <- broom::tidy(mod)
          # grab only the col2 row
          mod.sum <- mod.sum[mod.sum$term==col2,]

      } else if(reg.method == "lmer"){
          # Some data will have issues, this catches that error and returns NA
          tryCatch({
             mod <- run_single_lmer(paste0(col1, "~", covariates, " + ", col2),
                                   data)
             mod.sum <- broom.mixed::tidy(mod)

          }, error=function(e) {
              if(e$message == "not a positive definite matrix (and positive semidefiniteness is not checked)"){
                  mod.sum <- data.frame("effect"=NA, "group"=NA,
                                         "term"=col2, "estimate"=NA,
                                          "std.error"=NA, "statistic"=NA,
                                          "df"=NA, "p.value"=NA,
                                          "feature"=col1)
              }
          }) # END trycatch block
          # grab only the col2 beta row
          mod.sum <- mod.sum[mod.sum$term==col2 & mod.sum$effect=="fixed",]
      }

      mod.sum$feature <- col1

      mod.sum
  }

    if(show_progress){ close(pb) } # close the progress bar
    parallel::stopCluster(cl) # close the cluster for parallel computation



    mod.summaries <- as.data.frame(mod.summaries)
    # adjust p value only for non-intercept terms
    # and if we have a joint p, only adjust it for the beta coefficient (it's copied for the logistic)
    if((reg.method == "zibr") & zero_prop_from_formula){
        mod.summaries[mod.summaries$term!="(Intercept)" & mod.summaries$parameter == "beta",
                      "q"] <- p.adjust(
                          mod.summaries[mod.summaries$term!="(Intercept)" & mod.summaries$parameter == "beta",
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
        mod.summaries[mod.summaries$term!="(Intercept)" & mod.summaries$effect == "fixed",
                      "q"] <- p.adjust(
                          mod.summaries[mod.summaries$term!="(Intercept)"  & mod.summaries$effect == "fixed",
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
        stop(paste0(colnames.intersection, " is found in both datasets."))
    }
}


