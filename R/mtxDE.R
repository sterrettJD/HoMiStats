library(gamlss)

#' Run a single zero-inflatd beta regression using gamlss
#' @description runs a zero-inflated beta distribution GAM using gamlss
#' @param formula the full formula to be used in the regression
#' @param data A dataframe to be used in the regression
#' @return the resulting model
#' @export
#' @importFrom gamlss.dist BEZI
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
#'
run_single_beta_reg_zibr <- function(logistic_cov=NULL, beta_cov, Y,
                                     subject_ind=NULL, time_ind=NULL, data){
    # If no logistic covariate is passed, make constant data for fitting just an intercept
    if(is.null(logistic_cov)){
        logistic_cov_dat <- as.matrix(rep(0, nrow(data)))
    } else {
        logistic_cov_dat <- as.matrix(data[,logistic_cov])
    }

    # If no subject ID, make them up, assuming each sample comes from a unique participant
    if(is.null(subject_ind)){
        subject_ind_dat <- as.matrix(1:nrow(data))
    } else {
        subject_ind_dat <- as.matrix(data[,subject_ind])
    }
    # If no timepoint column, assume each sample is from the same timepoint
    if(is.null(time_ind)){
        time_ind_dat <- as.matrix(rep(1, nrow(data)))
    } else {
        time_ind_dat <- as.matrix(data[,time_ind])
    }


    mod <- ZIBR::zibr(
               logistic_cov=logistic_cov_dat,
               beta_cov=as.matrix(data[,beta_cov]),
               Y=as.matrix(data[,Y]),
               subject_ind=subject_ind_dat,
               time_ind=time_ind_dat,
               verbose=T)

    return(mod)
}

#' Tidy ZIBR output
#' @description ZIBR returns a list with lots of different kinds of entrie
#' @param mod A ZIBR model results object.
#' @return a dataframe with the tidy model results
#' @export
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
        print(term)
        clean[clean$term==term, "joint.p"] <- mod$joint_p[term]
    }

    return(clean)

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


#' Get random effects
#' @description Sometimes you want the name of random effects variables from a formula. This function does that.
#' @param form A formula with random effects.
#' @return a character vector with the random effects variable names, or NULL if there are no random effects.
#' @importFrom lme4 findbars
#' @export
#'
get_random_fx <- function(form){
    vars <- sapply(lme4::findbars(form),
                    FUN=function(x) as.character(x)[3])

    if(length(vars) == 0){
        return(NULL)
    }

    return(vars)

}



#' Run differential expression analysis for metatranscriptomics data
#' @description Regresses each feature using the formula provided and returns the statistical summary for each feature
#' @param formula the right hand side of the regression formula to be used for each feature
#' @param feature.table A dataframe, where rows are samples, and columns are genes/features. Row names should be sample IDs.
#' @param metadata The metadata dataframe, where rows are samples and columns are features
#' @param sampleID A string denoting name of the column in metadata with sample IDs, which should be used to merge metadata with the feature table's rownames
#' @param reg.method A string denoting the method to use for regression. Options include "zibr" and "gamlss".
#' @param padj A string denoting the p value adustment method. Options can be checked using 'p.adjust.methods'
#' @param zero_prop_from_formula In ZIBR zero-inflated beta regression, should the zeroes be modeled with the provided formula? Default is TRUE.
#' @param zibr_zibr_time_ind A string denoting the name of the time column for ZIBR. Defaults to NULL, which is implemented as a constant time value in ZIBR to not fit a time effect. This argument does nothing if reg.method is not "zibr".
#' @return a dataframe with differential expression results
#' @export
#' @importFrom broom.mixed tidy
#' @importFrom lme4 findbars
#'
run_mtxDE <- function(formula, feature.table, metadata, sampleID,
                      reg.method="zibr", padj="fdr",
                      zero_prop_from_formula=T,
                      zibr_time_ind=NULL){
    # Check for values of one, which the beta regression can't handle
    check_for_ones(feature.table)

    # merge the feature table and metadata based on the rownames
    formula <- paste0(" ~ ", formula)
    metadata.vars <- c(all.vars(as.formula(formula)), sampleID)
    data <- merge(feature.table, metadata[,metadata.vars], by="row.names")

    # Handle random effects
    # Extracts random effects from a formula
    random.effects.vars <- get_random_fx(as.formula(formula))

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
    }

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
        if(reg.method == "gamlss"){
            mod <- run_single_beta_reg_gamlss(paste0(col, formula),
                                              data=data)
            mod.sum <- broom.mixed::tidy(mod)
        }

        if((reg.method == "zibr") & zero_prop_from_formula){
            vars <- all.vars(as.formula(formula))
            mod <- run_single_beta_reg_zibr(logistic_cov=vars, beta_cov=vars,
                                            Y=col,
                                            subject_ind=NULL, time_ind=zibr_time_ind,
                                            data=data)

            mod.sum <- tidy_zibr_results(mod)
            }

        if((reg.method == "zibr") & zero_prop_from_formula==FALSE){
            vars <- all.vars(as.formula(formula))
            mod <- run_single_beta_reg_zibr(logistic_cov=NULL, beta_cov=vars,
                                            Y=col,
                                            subject_ind=NULL, time_ind=zibr_time_ind,
                                            data=data)

            mod.sum <- tidy_zibr_results(mod)
        }

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
