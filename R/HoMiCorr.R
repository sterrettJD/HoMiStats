
#' Create a correlation network between host and microbial data
#' @description Runs regressions between all features and returns a correlation matrix
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
#' @return a correlation matrix
#' @export
#' @importFrom broom.mixed tidy
#' @importFrom broom tidy
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel stopCluster
#' @importFrom doSNOW registerDoSNOW
#'
run_mtxDE <- function(mtx, host,
                      covariates=NULL, metadata=NULL,
                      sampleID=NULL,
                      reg.method="zibr",
                      padj="fdr",
                      zero_prop_from_formula=T,
                      zibr_time_ind=NULL,
                      ncores=1,
                      show_progress=TRUE
){
    # Check for values of one, which the beta regression can't handle
    if(reg.method %in% c("zibr", "gamlss")){
        check_for_ones(mtx)
    }

    # merge the feature table and metadata based on the rownames
    formula <- paste0(" ~ ", covariates)
    metadata.vars <- c(all.vars(as.formula(formula)), sampleID)

    data <- merge(mtx, host,
                  metadata[,c(metadata.vars, zibr_time_ind)],
                  by.x="row.names", by.y=sampleID)

    # initialize the model summaries df
    if(reg.method == "gamlss"){
        mod.summaries <- data.frame(matrix(nrow=0, ncol=7))
        colnames(mod.summaries) <- c("parameter", "term",
                                     "estimate", "std.error",
                                     "statistic", "p.value",
                                     "feature")

    } else if(reg.method == "zibr"){
        mod.summaries <- data.frame(matrix(nrow=0, ncol=6))
        colnames(mod.summaries) <- c("parameter", "term",
                                     "estimate",
                                     "p.value", "joint.p",
                                     "feature")
        # extract vars for zibr
        vars <- all.vars(as.formula(formula))
        # Extracts random effects from formula
        random.effects.vars <- get_random_fx(as.formula(formula))
        fixed.vars <- setdiff(vars, random.effects.vars)

    } else if(reg.method=="lm"){
        mod.summaries <- data.frame(matrix(nrow=0, ncol=6))
        colnames(mod.summaries) <- c("term", "estimate",
                                     "std.error", "statistic",
                                     "p.value", "feature")

    } else if(reg.method=="lmer"){
        mod.summaries <- data.frame(matrix(nrow=0, ncol=9))
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
        # Initializes progress bar
        total.features <- ncol(mtx)*ncol(host)
        n_iter <- (total.features*(total.features-1))/2
        i <- 0

        pb <- txtProgressBar(max=n_iter-1, style=3)
        progress <- function(n) setTxtProgressBar(pb, n)
        doSNOWopts <- list(progress = progress)
    } else {
        doSNOWopts <- list()
    }

    all.featurenames <- c(colnames(mtx), colnames(host))
    featurenames.combos <- combn(all.featurenames, 2,
                                 simplify=FALSE)

    # Loop through each combination of columns and run the regression
    mod.summaries <- foreach::foreach(cols=featurenames.combos,
                                      .combine=rbind,
                                      .options.snow = doSNOWopts) %dopar% {
        col1 <- cols[1]
        col2 <- cols[2]

      if(reg.method == "gamlss"){
          mod <- run_single_beta_reg_gamlss(paste0(col1, formula, " + ", col2),
                                            data=data)
          mod.sum <- broom.mixed::tidy(mod)
      } else if((reg.method == "zibr") & (zero_prop_from_formula==TRUE)){
          mod <- run_single_beta_reg_zibr(logistic_cov=fixed.vars, beta_cov=fixed.vars,
                                          Y=col,
                                          subject_ind=random.effects.vars,
                                          time_ind=zibr_time_ind,
                                          data=data)

          mod.sum <- tidy_zibr_results(mod)

      } else if((reg.method == "zibr") & (zero_prop_from_formula==FALSE)){

          mod <- run_single_beta_reg_zibr(logistic_cov=NULL, beta_cov=c(fixed.vars, col2),
                                          Y=col1,
                                          subject_ind=random.effects.vars,
                                          time_ind=zibr_time_ind,
                                          data=data)

          mod.sum <- tidy_zibr_results(mod)

      } else if(reg.method == "lm"){
          mod <- run_single_lm(paste0(col1, formula, " + ", col2), data)
          mod.sum <- broom::tidy(mod)

      } else if(reg.method == "lmer"){
          mod <- run_single_lmer(paste0(col, formula, " + ", col2), data)
          mod.sum <- broom.mixed::tidy(mod)
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
