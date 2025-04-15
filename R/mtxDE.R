#' Run a single zero-inflated beta regression using gamlss
#' @description runs a zero-inflated beta distribution GAM using gamlss
#' @param formula the full formula to be used in the regression
#' @param data A dataframe to be used in the regression
#' @param controller The `gamlss::gamlss.control()` object to be used by
#' the regression
#' @return the resulting model
#' @export
#' @importFrom gamlss.dist BEZI
#' @importFrom gamlss gamlss
#' @importFrom gamlss gamlss.control
#' @examples
#' data <- data.frame(a=c(0.1, 0.2, 0.3, 0.4),
#'                    b=c(0, 0, 1, 1))
#'
#' run_single_beta_reg_gamlss("a ~ b",
#'                            data=data)
#'
run_single_beta_reg_gamlss <- function(formula,
                                        data,
                                        controller=gamlss::gamlss.control(
                                                trace=FALSE)){
    mod <- gamlss::gamlss(as.formula(formula),
                        data=data,
                        family=gamlss.dist::BEZI,
                        control=controller)
    return(mod)
}


#' Run a single zero-inflated beta regression using gamlss
#' @description runs a zero-inflated beta regression using ZIBR
#' @param logistic_cov a string (or vector of strings)
#' indicating the variables to be used as covariates
#' predicting whether a value is 0 (in the logistic portion of ZIBR).
#' Defaults to NULL,
#' for which just an intercept will be fit for the logistic portion of ZIBR.
#' @param beta_cov a string (or vector of strings)
#' indicating the variables to be used as covariates predicting
#' the relative abundance of a non-zero value (in the beta portion of ZIBR)
#' @param Y a string indicating the the column
#' with zero-inflated relative abundance data
#' @param subject_ind a string indicating the column
#' with participant ID information for longitudinal data.
#' Defaults to NULL,
#' which creates unique participant ID for each sample
#' in the case of non-longitudinal, cross-sectional data
#' @param time_ind a string indication the column with timepoint data.
#' Defaults to NULL, which creates a constant data column for non-longitudinal,
#' cross-sectional data
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
    # If no logistic covariate is passed,
    # make constant data for fitting just an intercept
    if(is.null(logistic_cov)){
        logistic_cov_dat <- as.matrix(rep(0, nrow(data)))
    } else {
        logistic_cov_dat <- as.matrix(data[,logistic_cov])
        colnames(logistic_cov_dat) <- logistic_cov
    }

    # If no subject ID, make them up,
    # assuming each sample comes from a unique participant
    if(is.null(subject_ind)){
        subject_ind_dat <- as.matrix(seq_len(nrow(data)))
    } else {
        subject_ind_dat <- as.matrix(data[,subject_ind])
    }
    # If no timepoint column, assume each sample is from the same timepoint
    if(is.null(time_ind)){
        # If there are multiple timepoints for any participant
        if(length(unique(subject_ind_dat)) != length(subject_ind_dat)){
            stop("A timepoint column is necessary ",
            "if there are longitudinal samples.")
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
#' @description ZIBR returns a list with lots of different kinds of entries.
#' This tidies that up.
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
#' @description ZIBR outputs terms as "var1", "var2", etc in the results.
#' This maps these to the original variable names
#' @param data a character vector of the ZIBR terms
#' (e.g. c("intercept", "var1", "var2"))
#' @param real.names names of variables in the order they were passed to ZIBR.
#' Do not include intercept
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
#' @description The zero-inflated beta regression can't handle values of 1.
#' This checks for them and raises an error if they exist.
#' @param feature.table A dataframe, where rows are samples,
#' and columns are genes/features. Row names should be sample IDs.
#' @return Nothing
#' @export
#' @examples
#' data <- data.frame(a=c(0.1, 0.2, 0.3, 0.4),
#'                  b=c(0, 0, 1, 1))
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
#' @description The beta regression methods don't deal well
#' with features that are entirely undetected.
#' This filters them from the feature table.
#' @param feature.table A dataframe, where rows are samples,
#' and columns are genes/features. Row names should be sample IDs.
#' @return A dataframe, with undetected features removed
#' (where rows are samples, and columns are genes/features).
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
#' @description Gets the names of random effects variables from a formula.
#' @param form A formula with random effects.
#' @return a character vector with the random effects variable names,
#' or NULL if there are no random effects.
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

#' Check data validity for mtxDE function
#' @description This function checks the validity of your
#' data for mtxDE analyses. This will check that your
#' dna and feature tables have row names that match the given
#' sampleID metadata column. It also ensures that dna.table and
#' feature.table share column names.
#' @param feature.table A data frame where rows are samples and
#' columns are features (e.g., genes).
#' Row names should correspond to sample IDs.
#' @param dna.table A data frame where rows are samples and
#' columns are features (e.g., genes).
#' Row names should correspond to sample IDs.
#' This table should contain gene abundance data.
#' @param metadata A data frame containing metadata for the samples,
#' where rows are samples
#' and columns are metadata variables (e.g., phenotype, timepoint).
#' @param sampleID A string representing the column name in `metadata`
#' that contains the sample IDs.
#' This column is used to merge the metadata with the feature table.
#' @return Nothing.
#' @export
#' @examples
#' feature.table <- data.frame(a=c(0.1, 0.0, 0.0, 0.0),
#'                              b=c(0.5, 0.5, 0.5, 0.4),
#'                              c=c(0.4, 0.5, 0.0, 0.0),
#'                              d=c(0.0, 0.0, 0.5, 0.6))
#' dna.table <- data.frame(a=c(0.3, 0.1, 0.0, 0.0),
#'                              b=c(0.2, 0.2, 0.2, 0.1),
#'                              c=c(0.4, 0.5, 0.1, 0.0),
#'                              d=c(0.0, 0.0, 0.5, 0.6))
#' row.names(feature.table) <- paste0("sample_", seq_len(4))
#' row.names(dna.table) <- paste0("sample_", seq_len(4))
#' metadata <- data.frame(SampleID=paste0("sample_", seq_len(4)),
#'                           phenotype=c(0,0,1,1),
#'                           participant=c(0,1,0,1),
#'                           timepoint=c(0,0,1,1))
#' # This will not raise an error
#' check_data_mtxDE(feature.table = feature.table,
#'                  dna.table = dna.table,
#'                  metadata = metadata,
#'                  sampleID="SampleID")
#' # This will raise an error
#' feature.table <- data.frame(e=c(0.1, 0.0, 0.0, 0.0),
#'                             f=c(0.5, 0.5, 0.5, 0.4),
#'                             g=c(0.4, 0.5, 0.0, 0.0),
#'                             h=c(0.0, 0.0, 0.5, 0.6))
#' \dontrun{
#' check_data_mtxDE(feature.table = feature.table,
#'                  dna.table = dna.table,
#'                  metadata = metadata,
#'                  sampleID="SampleID")
#'}
check_data_mtxDE <- function(feature.table, dna.table=NULL, metadata, sampleID){

  if ((!is.null(metadata)) && (!is.null(sampleID))) {
    if ((sampleID %in% colnames(metadata)) == FALSE) {
      stop("sampleID column '", sampleID, "' not found in metadata.")
    }
    if (length(setdiff(rownames(feature.table),
                       metadata[, sampleID])) > 0) {
      stop("Row names of `feature.table` do not match `metadata$",
           sampleID, "`.")
    }
    if (!is.null(dna.table)){
      if (length(setdiff(rownames(dna.table),
                         metadata[, sampleID])) > 0) {
        stop("Row names of `dna.table` do not match `metadata$",
             sampleID, "`.")
      }
      if (length(intersect(colnames(dna.table),
                           colnames(feature.table))) < 1) {
        stop("None of the colnames of `dna.table` match `feature.table`")
      }
    }
  }
}

#' Filter two tables by shared columns
#' @description This function filters two data frames
#' to only retain shared columns.
#' @param table1 A data frame to be filtered.
#' @param table2 A second data frame to be filtered.
#' @param table1_name A character string specifying the name of `table1`.
#' @param table2_name A character string specifying the name of `table2`
#' @return A list containing:
#' \describe{
#'   \item{`table1_name`}{The filtered version of `table1`,
#'   keeping only the shared columns.}
#'   \item{`table2_name`}{The filtered version of `table2`,
#'   keeping only the shared columns.}
#' }
#' @export
#' @examples
#' feature.table <- data.frame(a = c(0.1, 0.2), b = c(0.3, 0.4), c = c(0.5, 0.6))
#' dna.table <- data.frame(a = c(0.7, 0.8), b = c(0.9, 1.0))
#' result <- filter_tables_by_shared_columns(feature.table, dna.table,
#'                                           "feature.table", "dna.table")
#' filtered_feature <- result$feature.table
#' filtered_dna <- result$dna.table
#'
filter_tables_by_shared_columns <- function(table1, table2, table1_name,
                                            table2_name) {
  all.feature.vars <- intersect(colnames(table1), colnames(table2))

  if (length(all.feature.vars) == 0) {
    stop("No shared columns between `feature.table` and `dna.table`.
         Ensure the tables have overlapping features.")
  }

  removed_from_table1 <- setdiff(colnames(table1), all.feature.vars)
  removed_from_table2 <- setdiff(colnames(table2), all.feature.vars)

  if (length(removed_from_table1) > 0) {
    warning(sprintf("The following feature(s) were removed from `%s`
                    because they do not have matching column(s) in `%s`: %s",
                    table1_name, table2_name,
                    paste(removed_from_table1, collapse = ", ")))
  }

  if (length(removed_from_table2) > 0) {
    warning(sprintf("The following feature(s) were removed from `%s`
                    because they do not have matching column(s) in `%s`: %s",
                    table2_name, table1_name,
                    paste(removed_from_table2, collapse = ", ")))
  }

  table1 <- table1[, all.feature.vars, drop = FALSE]  # Filter table1
  table2 <- table2[, all.feature.vars, drop = FALSE]  # Filter table2

  return(stats::setNames(list(table1, table2), c(table1_name, table2_name)))
}

#' Prepare Data for Differential Expression Analysis (internal)
#'
#' @description This function prepares the data by merging the feature table
#' with metadata based on the sample IDs.
#' It ensures that the sample IDs in the feature table match the
#' sample IDs in the metadata and validates the input for any inconsistencies.
#'
#' @param feature.table A data frame where rows are samples and
#' columns are features (e.g., genes).
#' Row names should correspond to sample IDs.
#' @param metadata A data frame containing metadata for the samples,
#' where rows are samples
#' and columns are metadata variables (e.g., phenotype, timepoint).
#' @param formula A string representing the right-hand side of the
#' regression formula to be used
#' in the analysis. All variables in this formula should be numeric.
#' @param sampleID A string representing the column name in `metadata`
#' that contains the sample IDs.
#' This column is used to merge the metadata with the feature table.
#' @param zibr_time_ind A string representing the time column in `metadata` for
#' Zero-Inflated Beta Regression (ZIBR).
#' This is used only when the regression method involves ZIBR.
#' @param dna.table A data frame where rows are samples and
#' columns are features (e.g., genes).
#' Row names should correspond to sample IDs.
#' This table should contain gene abundance data.
#'
#' @return A merged data frame containing the feature table and metadata,
#' ready for regression analysis.
#'
#' @keywords internal
#'
.prepare_data_mtxDE <- function(feature.table, metadata, formula,
                                sampleID, zibr_time_ind, dna.table=NULL) {
  # keep metadata columns that are in the formula
  metadata.vars <- c(all.vars(as.formula(formula)), sampleID)
  # check the validity of the data
  check_data_mtxDE(feature.table, dna.table, metadata, sampleID)
  if (!is.null(dna.table)) {
    # filter data tables to contain only shared columns
    filt.tables <- filter_tables_by_shared_columns(dna.table, feature.table,
                                                   "dna.table", "feature.table")
    dna.table <- filt.tables$dna.table
    feature.table <- filt.tables$feature.table
    # update dna.table col names before merge
    colnames(dna.table) <- paste0(colnames(dna.table), "_mgx")
    # merge feature table, dna table, and metadata based on row names
    data <- merge(feature.table,
                  metadata[,c(metadata.vars, zibr_time_ind)],
                  by.x="row.names", by.y=sampleID)
    data <- merge(data, dna.table, by.x="Row.names", by.y="row.names")
  } else {
    # merge feature table with metadata based on rownames
    data <- merge(feature.table,
                  metadata[,c(metadata.vars, zibr_time_ind)],
                  by.x="row.names", by.y=sampleID)
  }
  return(data)
}

#' Add matching dna column to formula or fixed vars
#'
#' @description
#' This function adds the matching dna feature of a `col` to a formula or
#' fixed vars
#'
#' @param data A data frame containing the merged feature table and metadata,
#' ready for regression.
#' @param col A string representing the name of the feature (column)
#' in `data` to be analyzed.
#' @param formula A string representing the formula for the regression model.
#' @param fixed.vars A vector of strings representing the fixed effect
#' variables to include in the regression.
#' @param reg.method A string indicating the regression method to be used.
#' Options include "zibr", "gamlss", "lm", and "lmer".
#'
#' @return list of updated `formula` and `fixed.vars` with the matching dna
#' column added to them
#'
#' @keywords internal
.add_dna_to_formula <- function(data, col, formula, fixed.vars, reg.method) {
  dna.col <- paste0(col, "_mgx")
  if (!(dna.col %in% colnames(data))) {
    stop("The following DNA feature was not found in metadata: ", dna.col)
  }
  if (reg.method == "zibr") {
    fixed.vars <- c(fixed.vars, dna.col)
  } else {
    formula <- paste0(formula, " + ", dna.col)
  }
  return(list(formula = formula, fixed.vars = fixed.vars))
}

#' Run a Single Regression for a Feature (internal)
#'
#' @description This function runs a single regression model for
#' a given feature using the specified regression method.
#' It handles different types of regression methods, including
#' Zero-Inflated Beta Regression (ZIBR),
#' Generalized Additive Models for Location Scale and Shape
#' (GAMLSS),linear regression (LM), and linear mixed effects regression (LME).
#'
#' @param data A data frame containing the merged feature table and metadata,
#' ready for regression.
#' @param reg.method A string indicating the regression method to be used.
#' Options include "zibr", "gamlss", "lm", and "lmer".
#' @param col A string representing the name of the feature (column)
#' in `data` to be analyzed.
#' @param formula A string representing the formula for the regression model.
#' @param fixed.vars A vector of strings representing the fixed effect
#' variables to include in the regression.
#' @param random.effects.vars A vector of strings representing the
#' random effect variables to include in the regression.
#' @param zibr_time_ind A string representing the time column in `data`
#' for ZIBR, or NULL if not applicable.
#' @param zero_prop_from_formula A boolean indicating whether zeroes should
#' be modeled using the provided formula (for ZIBR).
#' @param dna.table A data frame where rows are samples and
#' columns are features (e.g., genes).
#' Row names should correspond to sample IDs.
#' This table should contain gene abundance data.
#'
#' @return A data frame containing the summary of the regression results
#' for the specified feature.
#' The summary includes terms, estimates, standard errors, statistics,
#' p-values, and feature information.
#'
#' @keywords internal
#'
.run_single_regression_mtxDE <- function(data, reg.method,
                                        col, formula,
                                        fixed.vars=NULL,
                                        random.effects.vars=NULL,
                                        zibr_time_ind=NULL,
                                        zero_prop_from_formula=NULL,
                                        dna.table=NULL) {
   if (!is.null(dna.table)) {
    # check that dna col exists and add to formula/fixed.vars
    updated.dna <- .add_dna_to_formula(data, col, formula,
                                      fixed.vars, reg.method)
    formula <- updated.dna$formula
    fixed.vars <- updated.dna$fixed.vars
  }
    if (reg.method == "gamlss"){
        mod <- run_single_beta_reg_gamlss(paste0(col, formula),
                                            data=data)
        mod.sum <- broom.mixed::tidy(mod)
    } else if ((reg.method == "zibr") && (zero_prop_from_formula==TRUE)) {
        mod <- run_single_beta_reg_zibr(logistic_cov=fixed.vars,
                                        beta_cov=fixed.vars, Y=col,
                                        subject_ind=random.effects.vars,
                                        time_ind=zibr_time_ind, data=data)
        mod.sum <- tidy_zibr_results(mod)
        mod.sum$term <- map_zibr_termnames(mod.sum$term, fixed.vars)
    } else if ((reg.method == "zibr") && (zero_prop_from_formula==FALSE)) {
        mod <- run_single_beta_reg_zibr(logistic_cov=NULL,
                                        beta_cov=fixed.vars, Y=col,
                                        subject_ind=random.effects.vars,
                                        time_ind=zibr_time_ind, data=data)
        mod.sum <- tidy_zibr_results(mod)
        mod.sum$term <- map_zibr_termnames(mod.sum$term, fixed.vars)
    } else if (reg.method == "lm") {
        mod <- run_single_lm(paste0(col, formula), data)
        mod.sum <- broom::tidy(mod)
    } else if (reg.method == "lmer") {
        mod <- run_single_lmer(paste0(col, formula), data)
        mod.sum <- broom.mixed::tidy(mod)
    }
    mod.sum$feature <- col
    return(mod.sum)
}

#' Run Multiple Regressions for Features in Parallel (internal)
#'
#' @description This function runs regressions for each feature (column)
#' in the `feature.table` using the specified regression method.
#' The function handles different regression methods such as
#' Zero-Inflated Beta Regression (ZIBR),
#' Generalized Additive Models for Location Scale and Shape (GAMLSS),
#' linear regression (LM), and linear mixed effects regression (LME).
#'
#' @param data A data frame containing the merged feature table
#' and metadata, prepared for regression.
#' @param reg.method A string indicating the regression method to be used.
#' Options include
#' "zibr", "gamlss", "lm", and "lmer".
#' @param feature.table A data frame where rows represent samples and
#' columns represent features.
#' @param formula A string representing the formula for the regression model.
#' @param zero_prop_from_formula A boolean indicating whether zeroes
#' should be modeled using the provided formula (for ZIBR).
#' @param zibr_time_ind A string representing the time column in `data`
#' for ZIBR, or NULL if not applicable.
#' @param ncores An integer specifying the number of cores to use
#' for parallel processing.
#' @param show_progress A boolean indicating whether a progress bar
#' should be shown during parallel computation.
#' @param dna.table A data frame where rows are samples and
#' columns are features (e.g., genes).
#' Row names should correspond to sample IDs.
#' This table should contain gene abundance data.
#'
#' @return A data frame containing the regression summaries for all features.
#' The summaries include terms,
#' estimates, standard errors, statistics, p-values, and feature information.
#'
#' @keywords internal
#'
.run_regressions_mtxDE <- function(data, reg.method,
                                    feature.table, formula,
                                    zero_prop_from_formula,
                                    zibr_time_ind, ncores, show_progress,
                                   dna.table=NULL){
    if(reg.method=="zibr"){
        form.as.form <- as.formula(formula)
        random.effects.vars <- get_random_fx(form.as.form)
        fixed.vars <- setdiff(all.vars(form.as.form), random.effects.vars)
    } else {
      fixed.vars <- NULL
    }
    # Setup cluster for parallel running
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl, cores=ncores)
    doSNOW::registerDoSNOW(cl)
    if(show_progress){
        # Initializes progress bar
        n_iter <- ncol(feature.table)
        pb <- txtProgressBar(max=n_iter-1, style=3)
        progress <- function(n) setTxtProgressBar(pb, n)
        doSNOWopts <- list(progress = progress)
    } else {
        doSNOWopts <- list()
    }
    # Loop through each column and run the regression
    mod.summaries <- foreach::foreach(col=colnames(feature.table),
                                        .combine=rbind,
                                        .options.snow = doSNOWopts,
                                        .packages = "HoMiStats"
                                      ) %dopar% {
    .run_single_regression_mtxDE(data, reg.method,
                                col, formula,
                                fixed.vars, random.effects.vars,
                                zibr_time_ind,
                                zero_prop_from_formula,
                                dna.table=dna.table)
    }

    if(show_progress){ close(pb) } # close the progress bar
    parallel::stopCluster(cl) # close the cluster for parallel computation
    mod.summaries <- as.data.frame(mod.summaries)
    return(mod.summaries)
}


#' Run differential expression analysis for metatranscriptomics data
#' @description Regresses each feature using the formula provided and
#' returns the statistical summary for each feature
#' @param formula the right hand side of the regression formula to be used
#' for each feature.
#' All variables in this formula should be of a numeric type
#' (factors should be dummy coded, such that the reference value is 0).
#' @param feature.table A dataframe, where rows are samples,
#' and columns are genes/features.
#' Row names should be sample IDs.
#' @param metadata The metadata dataframe,
#' where rows are samples and columns are features
#' @param sampleID A string denoting name of the column in metadata
#' with sample IDs,
#' which should be used to merge metadata with the feature table's rownames
#' @param reg.method A string denoting the method to use for regression.
#' Options include "zibr" (Zero-inflated beta regression with random effects),
#' "gamlss" (Zero-inflated beta regression implemented via GAMLSS),
#' "lm" (linear regression), and
#' "lmer" (linear mixed effects regression implemented via lme4 and lmerTest).
#' @param padj A string denoting the p value adustment method.
#' Options can be checked using 'p.adjust.methods'
#' @param zero_prop_from_formula In ZIBR zero-inflated beta regression,
#' should the zeroes be modeled with the provided formula? Default is TRUE.
#' @param zibr_time_ind A string denoting the name of the time column
#' for ZIBR.
#' Defaults to NULL, which is implemented as a constant time value in ZIBR
#' to not fit a time effect.
#' This argument does nothing if reg.method is not "zibr".
#' @param ncores An integer denoting the number of cores to use
#' if running in parallel. Defaults to 1 (not parallelized).
#' @param show_progress A boolean denoting if a progress bar should be shown.
#' @return a dataframe with differential expression results
#' @param dna.table A data frame where rows are samples and
#' columns are features (e.g., genes).
#' Row names should correspond to sample IDs.
#' This table should contain gene abundance data.
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
#'                      metadata,
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
                        show_progress=TRUE,
                        dna.table=NULL){
    if(reg.method %in% c("zibr", "gamlss")){
        check_for_ones(feature.table) # Beta regression can't handle ones
        feature.table <- filter_undetected(feature.table) # or undetected feats
        if(!is.null(dna.table)){
          check_for_ones(dna.table)
          dna.table <- filter_undetected(dna.table)
        }
    }
    formula <- paste0(" ~ ", formula)
    data <- .prepare_data_mtxDE(feature.table, metadata,
                                formula, sampleID, zibr_time_ind,
                                dna.table = dna.table)

    # Loop through each column and run the regression
    mod.summaries <- .run_regressions_mtxDE(data, reg.method,
                                            feature.table, formula,
                                            zero_prop_from_formula,
                                            zibr_time_ind,
                                            ncores, show_progress,
                                            dna.table=dna.table)
    # Calling this from HoMiCorr
    adjusted.mod.summaries <- .adjust_p_values(mod.summaries, reg.method,
                                                zero_prop_from_formula, padj)
    return(adjusted.mod.summaries)
}
