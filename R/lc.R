#' Fit latent-class MPT models
#'
#' This function estimates one or more latent-class MPT models.
#' Starting from a one-class solution, it increments the number of latent classes
#' until a criterion is met.
#'
#' @param model_file Character. File path to an .eqn file
#' @param data_file Character. File path to a data file
#' @param classes,max_classes The number of classes for a single model or the
#' maximum number of classes. To estimate a single model, specify \code{classes};
#' to estimate multiple models with increasing number of classes, specify \code{max_classes}.
#' @param runs Integer. Number of optimization runs.
#' @param fisher_information The type of Fisher Information to be computed.
#'    Can be either \code{expected}, \code{montecarlo}, \code{observed}, or \code{none}.
#'    Defaults to \code{expected}. However, if expected Fisher Information cannot be
#'    computed, the Monte Carlo method is used.
#' @param montecarlo_samples The number of simulations to be used for Monte Carlo
#'   Fisher Information. Defaults to 1e5.
#' @param crit If multiple models with an increasing number of classes are to be estimated,
#'   use this criterion to stop estimating more complex models. Can be one of \code{"AIC"},
#'   \code{"BIC"}, or \code{"BIC2"}.
#' @param verbose Logical. Indicating whether an announcement is printed on the
#'   console when fitting of a new model starts.
#' @param keep_files Logical. Should temporary files be retained?
#'
#' @references
#'   Stahl, C., & Klauer, K.C. (2007). HMMTree: A computer program for hierarchical
#'   multinomial processing tree models. \emph{Behavior Research Methods}, \emph{39}, 267-273.
#'
#' @export

lc <- function(
  model_file
  , data_file
  , classes = NULL
  , max_classes = 20
  , runs = 20
  , fisher_information = "expected"
  , montecarlo_samples = 1e5
  , crit = "AIC"
  , verbose = FALSE
  , keep_files = FALSE
){

  # Check if running on Windows
  if(!Sys.info()[["sysname"]]=="Windows") {
    stop("Sorry, but model estimation only works on Windows machines. We know this must be disappointing.")
  }



  # Input validation ----
  if (!(crit %in% c("AIC","BIC", "BIC2"))) {
    stop("Parameter crit must be one of 'AIC', 'BIC', or 'BIC2'.")
  }



  # Create legacy eqn file
  eqn_file <- paste0("HMMTreeR-tmpfile-", basename(model_file))
  eqn_conv_out <- simplify_eqn(model_file = model_file, eqn_file = eqn_file)


  # Create HMMTree data file (tab-separated, first column is dataset name)
  dataset_name <- gsub(basename(data_file), pattern = ".dat|.csv|.txt| ", replacement = "")

  dat_file <- paste0("HMMTreeR-tmpfile-", dataset_name, ".dat")

  if(grepl(x = data_file, pattern = ".csv")) {
    tmp_dat <- utils::read.csv(file = data_file)
    if(any(grepl(x = tmp_dat[, 1], pattern = ";"))) {
      tmp_dat <- utils::read.csv2(file = data_file, check.names = FALSE)
      # print(tmp_dat)
    }
  } else {
    if(grepl(x = data_file, pattern = ".dat")) {
      tmp_dat <- utils::read.delim(file = data_file, check.names = FALSE)
    }
  }
  # print(tmp_dat)
  if(!all(eqn_conv_out$cat_order %in% colnames(tmp_dat))) {
    stop("Some categories defined in the .eqn file are missing from data.")
  }

  # try_conv <- suppressWarnings(as.integer(eqn_conv_out$cat_order))
  # print(eqn_conv_out$cat_order)
  # print(colnames(tmp_dat))
  # if(any(is.na(try_conv))) {
    tmp_dat <- tmp_dat[, eqn_conv_out$cat_order]
    tmp_dat <- cbind(
      data.frame(
        id = dataset_name
      )
      , tmp_dat
    )
  # }

  utils::write.table(
    x = tmp_dat
    , file = dat_file
    , sep = "\t"
    , row.names = FALSE
    , col.names = TRUE
    , quote = FALSE
  )

  nsubj <- nrow(tmp_dat)

  res <- list()

  # test fit: is (co)variance accounted for by single-class model?
  #pval <- pchisq(q=as.numeric(res[[nclass]]$S1), df=as.numeric(res[[nclass]]$S1_df), lower.tail=FALSE)
  #if ((length(pval)>0) || (pval < .05)){ # aggregate does not fit, estimate latent-class models

  if(!is.null(classes)) {
    # fit only one solution
    res <- model_object(
      hmmtreec(
        model_file = eqn_file
        , data_file = dat_file
        , nsubj = nsubj
        , nclass = classes
        , nruns = runs
        , fi = switch(fisher_information, "expected" = 3L, "montecarlo" = 2L, "observed" = 1L, "none" = 0L)
        , mc = montecarlo_samples
        , keep_files = keep_files
      )
    )
  } else {
    # Estimate more complex models until criterion is reached
    n_classes <- 1L
    repeat {
      if(verbose) {
        if(n_classes == 1L) {
          message("Fitting a latent-class MPT with one class (as a baseline).")
        } else {
          message("Fitting a latent-class MPT with ", n_classes, " classes.")
        }
      }

      tmp_hmmtreec <- try(
        hmmtreec(
          model_file = eqn_file
          , data_file = dat_file
          , nsubj = nsubj
          , nclass = n_classes
          , nruns = runs
          , fi = switch(fisher_information, "expected" = 3, "montecarlo" = 2, "observed" = 1, "none" = 0)
          , mc = montecarlo_samples
          , keep_files = keep_files
        )
      )

      if(class(tmp_hmmtreec) == "try-error") {
        break
      }

      tmp <- model_object(tmp_hmmtreec)

      # Stop estimating more complex models if failcodes occur
      if(tmp$description$failcode!=0) {
        break
      }

      # Stop if negative degrees of freedom occur (seems to be a byproduct of ill-conditioned FI)
      if(any(tmp$fit_statistics[c("df_M1", "df_M2", "df_S1", "df_S2")]<0)) {
        break
      }

      # Stop estimating more complex models if criterion (i.e., AIC or BIC) increased
      if(n_classes > 1) {
        crit_new <- tmp$fit_statistics[[crit]]
        crit_old <- res[[n_classes-1]]$fit_statistics[[crit]]

        if(isTRUE(crit_new > crit_old)) break
      }

      res[[n_classes]] <- tmp

      n_classes <- n_classes + 1
      if(n_classes > max_classes) {
        break
      }
    }
    class(res) <- c("lc_mpt_list", class(res))
  }

  # remove temporary files
  if(!keep_files) {
    file.remove(eqn_file, dat_file)
  }

  # return
  res
}
