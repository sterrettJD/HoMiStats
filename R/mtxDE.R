#' Run a single zero-inflated beta regression using gamlss
#' @description runs a zero-inflated beta distribution GAM using gamlss
#' @param formula the full formula to be used in the regression
#' @param data A dataframe to be used in the regression
#' @return the resulting model
#' @export
#' @importFrom gamlss.dist BEZI
#' @examples
#' data <- data.frame(a=c(0.1, 0.2, 0.3, 0.4),
#'                    b=c(0, 0, 1, 1))
#'
#' run_single_beta_reg_gamlss("a ~ b",
#'                            data=data)
#'
run_single_beta_reg_gamlss <- function(formula, data,
                                              controller=gamlss.control(trace=FALSE)){
    mod <- gamlss::gamlss(as.formula(formula),
                        data=data,
                        family=gamlss.dist::BEZI,
                        control=controller)
    return(mod)
}


#' Run a single zero-inflated beta regression using gamlss
#' @description runs a zero-inflated beta regression using ZIBR
#' @param logistic_cov a string (or vector of strings) indicating the variables to be used as covariates predicting whether a value is 0 (in the logistic portion of ZIBR). Defaults to NULL, for which just an intercept will be fit for the logistic portion of ZIBR.
#' @param beta_cov a string (or vector of strings) indicating the variables to be used as covariates predicting the relative abundance of a non-zero value (in the beta portion of ZIBR)
#' @param Y a string indicating the the column with zero-inflated relative abundance data
#' @param subject_ind a string indicating the column with participant ID information for longitudinal data. Defaults to NULL, which creates unique participant ID for each sample in the case of non-longitudinal, cross-sectional data
#' @param time_ind a string indication the column with timepoint data. Defaults to NULL, which creates a constant data column for non-longitudinal, cross-sectional data
#' @param data A dataframe to be used in the regression
#' @return the resulting model
#' @export
#' @importFrom ZIBR zibr
#' @examples
#' data <- data.frame(y=c(0.1, 0, 0.3, 0.4),
#'                    treated=c(0, 0, 0, 1),
#'                    participant=c(1, 2, 1, 2),
#'                    timepoint=c(0, 0, 1, 1))
#'
#' # This runs ZIBR with a constant amount of zero-inflation
#' # for all participants and random intercepts for participant
#' run_single_beta_reg_zibr(beta_cov="treated", Y="y",
#'                          subject_ind="participant",
#'                          time_ind="timepoint",
#'                          data=data)
#'
#' # This runs ZIBR where treatment impacts both whether a feature might be zero
#' # and the abundance of that feature (and random intercepts for participant)
#' run_single_beta_reg_zibr(logistic_cov="treated",
#'                          beta_cov="treated",
#'                          Y="y",
#'                          subject_ind="participant",
#'                          time_ind="timepoint",
#'                          data=data)
#'
#' # This can also be run without random effects
#' # by not specifying participant or timepoint
#' run_single_beta_reg_zibr(beta_cov="treated", Y="y",
#'                          data=data)
#'
run_single_beta_reg_zibr <- function(logistic_cov=NULL, beta_cov, Y,
                                     subject_ind=NULL, time_ind=NULL, data){
    # If no logistic covariate is passed, make constant data for fitting just an intercept
    if(is.null(logistic_cov)){
        logistic_cov_dat <- as.matrix(rep(0, nrow(data)))
    } else {
        logistic_cov_dat <- as.matrix(data[,logistic_cov])
        colnames(logistic_cov_dat) <- logistic_cov
    }

    # If no subject ID, make them up, assuming each sample comes from a unique participant
    if(is.null(subject_ind)){
        subject_ind_dat <- as.matrix(seq_len(nrow(data)))
    } else {
        subject_ind_dat <- as.matrix(data[,subject_ind])
    }
    # If no timepoint column, assume each sample is from the same timepoint
    if(is.null(time_ind)){
        # If there are multiple timepoints for any participant
        if(length(unique(subject_ind_dat)) != length(subject_ind_dat)){
            stop("A timepoint column is necessary if there are longitudinal samples.")
        }
        time_ind_dat <- as.matrix(rep(1, nrow(data)))
    } else {
        time_ind_dat <- as.matrix(data[,time_ind])
    }

    beta_cov_dat <- as.matrix(data[,beta_cov])
    colnames(beta_cov_dat) <- beta_cov

    mod <- ZIBR::zibr(
               logistic_cov=logistic_cov_dat,
               beta_cov=beta_cov_dat,
               Y=as.matrix(data[,Y]),
               subject_ind=subject_ind_dat,
               time_ind=time_ind_dat,
               verbose=TRUE)

    return(mod)
}

#' Tidy ZIBR output
#' @description ZIBR returns a list with lots of different kinds of entries. This tidies that up.
#' @param mod A ZIBR model results object.
#' @return a dataframe with the tidy model results
#' @export
#' @examples
#' data <- data.frame(y=c(0.1, 0, 0.3, 0.4),
#'                    treated=c(0, 0, 0, 1),
#'                    participant=c(1, 2, 1, 2),
#'                    timepoint=c(0, 0, 1, 1))
#'
#' # This runs ZIBR where treatment impacts both whether a feature might be zero
#' # and the abundance of that feature (and random intercepts for participant)
#' mod <- run_single_beta_reg_zibr(logistic_cov="treated",
#'                                 beta_cov="treated",
#'                                 Y="y",
#'                                 subject_ind="participant",
#'                                 time_ind="timepoint",
#'                                 data=data)
#' tidied.mod <- tidy_zibr_results(mod)
#'
tidy_zibr_results <- function(mod){
    # There is a joint_p if logistic_cov and beta_cov are the same
    clean <- data.frame(
        parameter=c(rep("logistic", nrow(mod$logistic_est_table)),
                    rep("beta", nrow(mod$beta_est_table))),
        term=c(row.names(mod$logistic_est_table),
               row.names(mod$beta_est_table)),
        estimate=c(mod$logistic_est_table[,"Estimate"],
                   mod$beta_est_table[,"Estimate"]),
        p.value=c(mod$logistic_est_table[,"Pvalue"],
                  mod$beta_est_table[,"Pvalue"]),
        joint.p=rep(NA,
                    nrow(mod$logistic_est_table) + nrow(mod$beta_est_table))
        )

    for(term in names(mod$joint_p)){
        clean[clean$term==term, "joint.p"] <- mod$joint_p[term]
    }

    return(clean)

}


#' Map ZIBR term names
#' @description ZIBR outputs terms as "var1", "var2", etc in the results. This maps these to the original variable names
#' @param data a character vector of the ZIBR terms (e.g. c("intercept", "var1", "var2"))
#' @param real.names names of variables in the order they were passed to ZIBR. Do not include intercept
#' @return a character vector with the mapped variable names
#' @importFrom plyr mapvalues
#' @export
#' @examples
#' data <- data.frame(y=c(0.1, 0, 0.3, 0.4),
#'                    treated=c(0, 0, 0, 1),
#'                    participant=c(1, 2, 1, 2),
#'                    timepoint=c(0, 0, 1, 1))
#'
#' # This runs ZIBR where treatment impacts both whether a feature might be zero
#' # and the abundance of that feature (and random intercepts for participant)
#' mod <- run_single_beta_reg_zibr(logistic_cov="treated",
#'                                 beta_cov="treated",
#'                                 Y="y",
#'                                 subject_ind="participant",
#'                                 time_ind="timepoint",
#'                                 data=data)
#' tidied.mod <- tidy_zibr_results(mod)
#' tidied.mod$term <- map_zibr_termnames(tidied.mod$term, c("treated"))
#'
map_zibr_termnames <- function(data, real.names){
    names(real.names) <- paste0("var", seq_len(length(real.names)))
    real.names <- c(real.names, "intersept"="intercept")
    return(plyr::mapvalues(data, from=names(real.names), to=real.names))
}


#' Run a single linear regression
#' @description runs a simple linear regression
#' @param formula the full formula to be used in the regression
#' @param data A dataframe to be used in the regression
#' @return the resulting model
#' @export
#' @importFrom stats lm
#' @examples
#' data <- data.frame(a=c(0.1, 0.2, 0.3, 0.4),
#'                    b=c(0, 0, 1, 1))
#'
#' run_single_lm("a ~ b",
#'               data=data)
#'
run_single_lm <- function(formula, data){
    mod <- lm(as.formula(formula),
              data=data)
    return(mod)
}


#' Run a single linear mixed effects regression
#' @description runs a linear mixed regression using lme4
#' @param formula the full formula to be used in the regression
#' @param data A dataframe to be used in the regression
#' @return the resulting model
#' @export
#' @importFrom lmerTest lmer
#' @examples
#' data <- data.frame(a=c(0.1, 0.2, 0.3, 0.4),
#'                    b=c(0, 0, 1, 1),
#'                    participant=c(0, 1, 0, 1))
#'
#' run_single_lmer("a ~ b + (1|participant)",
#'               data=data)
#'
run_single_lmer <- function(formula, data){
    mod <- lmerTest::lmer(as.formula(formula),
              data=data)
    return(mod)
}


#' Check for 1 in feature.table
#' @description The zero-inflated beta regression can't handle values of 1. This checks for them and raises an error if they exist.
#' @param feature.table A dataframe, where rows are samples, and columns are genes/features. Row names should be sample IDs.
#' @return Nothing
#' @export
#' @examples
#' data <- data.frame(a=c(0.1, 0.2, 0.3, 0.4),
#'                    b=c(0, 0, 1, 1))
#'
#' # This would raise an error about rows 3 and 4
#' try(check_for_ones(data))
#'
check_for_ones <- function(feature.table){
    rows.with.ones <- which(feature.table >= 1, arr.ind = TRUE)[,1]
    unique.rows.with.ones <- unique(rows.with.ones)
    unique.rownames.with.ones <- rownames(feature.table)[unique.rows.with.ones]
    if(length(unique.rownames.with.ones) > 0){
        rownames.string <- paste(unique.rownames.with.ones, collapse=", ")
        stop("The following rows contains a value of one ",
             "which the zero-inflated beta regression cannot handle: ",
              rownames.string)
    }
}

#' Filter undetected features from the feature.table
#' @description The beta regression methods don't deal well with features that are entirely undetected. This filters them from the feature table.
#' @param feature.table A dataframe, where rows are samples, and columns are genes/features. Row names should be sample IDs.
#' @return A dataframe, with undetected features removed (where rows are samples, and columns are genes/features).
#' @export
#' @examples
#' good.data <- data.frame(a=c(0.1, 0.2, 0.3, 0.4),
#'                    b=c(0, 0, 1, 1))
#' bad.data <- data.frame(a=c(0.1, 0.2, 0.3, 0.4),
#'                    b=c(0, 0, 1, 1),
#'                    undetected_feature=c(0, 0, 0, 0))
#'
#' # This will not raise any warning
#' filter_undetected(good.data)
#'
#' # This would raise a warning
#' try(filter_undetected(bad.data))
#'
filter_undetected <- function(feature.table){
    detected.features <- which(colSums(feature.table) > 0)
    if(length(detected.features) < ncol(feature.table)){
        warning(strwrap(prefix = " ", initial = "",
                        "Undetected features (summing to 0)
                        were removed from this dataset.
                        Consider filtering these prior to
                        analysis."))
    }
    return(feature.table[, detected.features])
}


#' Get random effects
#' @description Sometimes you want the name of random effects variables from a formula. This function does that.
#' @param form A formula with random effects.
#' @return a character vector with the random effects variable names, or NULL if there are no random effects.
#' @importFrom lme4 findbars
#' @export
#' @examples
#' form <- as.formula("y ~ x + (1|z)")
#'
#' # This would return `c("z")`
#' get_random_fx(form)
#'
#' form.norand <- as.formula("y ~ x")
#' # This would return NULL
#' get_random_fx(form.norand)
#'
get_random_fx <- function(form){
    vars <- vapply(lme4::findbars(form),
                   FUN = function(x) as.character(x)[3],
                   FUN.VALUE = character(1))


    if(length(vars) == 0){
        return(NULL)
    }

    return(vars)
}


#' Run differential expression analysis for metatranscriptomics data
#' @description Regresses each feature using the formula provided and returns the statistical summary for each feature
#' @param formula the right hand side of the regression formula to be used for each feature. All variables in this formula should be of a numeric type (factors should be dummy coded, such that the reference value is 0).
#' @param feature.table A dataframe, where rows are samples, and columns are genes/features. Row names should be sample IDs.
#' @param metadata The metadata dataframe, where rows are samples and columns are features
#' @param sampleID A string denoting name of the column in metadata with sample IDs, which should be used to merge metadata with the feature table's rownames
#' @param reg.method A string denoting the method to use for regression. Options include "zibr" (Zero-inflated beta regression with random effects), "gamlss" (Zero-inflated beta regression implemented via GAMLSS), "lm" (linear regression), and "lmer" (linear mixed effects regression implemented via lme4 and lmerTest).
#' @param padj A string denoting the p value adustment method. Options can be checked using 'p.adjust.methods'
#' @param zero_prop_from_formula In ZIBR zero-inflated beta regression, should the zeroes be modeled with the provided formula? Default is TRUE.
#' @param zibr_zibr_time_ind A string denoting the name of the time column for ZIBR. Defaults to NULL, which is implemented as a constant time value in ZIBR to not fit a time effect. This argument does nothing if reg.method is not "zibr".
#' @param ncores An integer denoting the number of cores to use if running in parallel. Defaults to 1 (not parallelized).
#' @param show_progress A boolean denoting if a progress bar should be shown.
#' @return a dataframe with differential expression results
#' @export
#' @importFrom broom.mixed tidy
#' @importFrom broom tidy
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @examples
#' feature.table <- data.frame(a=c(0.1, 0.0, 0.0, 0.0),
#'                                b=c(0.5, 0.5, 0.5, 0.4),
#'                                c=c(0.4, 0.5, 0.0, 0.0),
#'                                d=c(0.0, 0.0, 0.5, 0.6))
#' row.names(feature.table) <- paste0("sample_", seq_len(4))
#' metadata <- data.frame(SampleID=paste0("sample_", seq_len(4)),
#'                           phenotype=c(0,0,1,1),
#'                           participant=c(0,1,0,1),
#'                           timepoint=c(0,0,1,1))
#'
#' zibr.no.random.fx <- run_mtxDE("phenotype",
#'                                feature.table,
#'                                metadata,
#'                                sampleID="SampleID",
#'                                show_progress=FALSE)
#'
#' # Use ZIBR with random effects
#' zibr.with.random.fx <- run_mtxDE("phenotype + (1|participant)",
#'                                  feature.table,
#'                                  metadata, sampleID="SampleID",
#'                                  zibr_time_ind="timepoint",
#'                                  show_progress=FALSE)
#'
#' # Use gamlss zero inflated beta regression with no random effects
#' gamlss.no.random.fx <- run_mtxDE("phenotype", feature.table,
#'                                   metadata, sampleID="SampleID",
#'                                   reg.method="gamlss",
#'                                   show_progress=FALSE)
#'
#' # Use a linear model
#' lm.res <- run_mtxDE("phenotype",
#'                      feature.table,
#                       metadata,
#'                      reg.method="lm",
#'                      sampleID="SampleID",
#'                      show_progress=FALSE)
#'
#' # Use a linear model with random effects for participant
#' lmer.res <- run_mtxDE("phenotype + (1|participant)",
#'                        feature.table,
#'                        metadata,
#'                        reg.method="lmer",
#'                        sampleID="SampleID",
#'                        show_progress=FALSE)
#'
run_mtxDE <- function(formula, feature.table, metadata, sampleID,
                      reg.method="zibr", padj="fdr",
                      zero_prop_from_formula=TRUE,
                      zibr_time_ind=NULL,
                      ncores=1,
                      show_progress=TRUE
                      ){
    # Check for values of one, which the beta regression can't handle
    if(reg.method %in% c("zibr", "gamlss")){
        check_for_ones(feature.table)
        # gamlss can't handle undetected features
        feature.table <- filter_undetected(feature.table)
    }


    # merge the feature table and metadata based on the rownames
    formula <- paste0(" ~ ", formula)
    metadata.vars <- c(all.vars(as.formula(formula)), sampleID)

    data <- merge(feature.table,
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
        n_iter <- ncol(feature.table)
        i <- 0

        pb <- txtProgressBar(max=n_iter-1, style=3)
        progress <- function(n) setTxtProgressBar(pb, n)
        doSNOWopts <- list(progress = progress)
    } else {
        doSNOWopts <- list()
    }

    # Loop through each column and run the regression
    mod.summaries <- foreach::foreach(col=colnames(feature.table),
                                      .combine=rbind,
                                      .options.snow = doSNOWopts) %dopar% {

        if(reg.method == "gamlss"){
            mod <- run_single_beta_reg_gamlss(paste0(col, formula),
                                              data=data)
            mod.sum <- broom.mixed::tidy(mod)
        } else if((reg.method == "zibr") & (zero_prop_from_formula==TRUE)){
            mod <- run_single_beta_reg_zibr(logistic_cov=fixed.vars, beta_cov=fixed.vars,
                                            Y=col,
                                            subject_ind=random.effects.vars,
                                            time_ind=zibr_time_ind,
                                            data=data)

            mod.sum <- tidy_zibr_results(mod)
            mod.sum$term <- map_zibr_termnames(mod.sum$term, fixed.vars)

        } else if((reg.method == "zibr") & (zero_prop_from_formula==FALSE)){

            mod <- run_single_beta_reg_zibr(logistic_cov=NULL, beta_cov=fixed.vars,
                                            Y=col,
                                            subject_ind=random.effects.vars,
                                            time_ind=zibr_time_ind,
                                            data=data)

            mod.sum <- tidy_zibr_results(mod)
            mod.sum$term <- map_zibr_termnames(mod.sum$term, fixed.vars)

        } else if(reg.method == "lm"){
            mod <- run_single_lm(paste0(col, formula), data)
            mod.sum <- broom::tidy(mod)

        } else if(reg.method == "lmer"){
            mod <- run_single_lmer(paste0(col, formula), data)
            mod.sum <- broom.mixed::tidy(mod)
        }

        mod.sum$feature <- col

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
