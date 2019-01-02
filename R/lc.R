#' Fit latent-class MPT models
#'
#' This function estimates one or more latent-class MPT models.
#' Starting from a one-class solution, it increments the number of latent classes
#' until a criterion is met.
#'
#' @param classes,max_classes The number of classes for a single model or the
#' maximum number of classes. To estimate a single model, specify \code{classes};
#' to estimate multiple models with increasing number of classes, specify \code{max_classes}.
#' @param fisher_information The type of Fisher Information to be computed.
#'    Can be either \code{expected}, \code{montecarlo}, \code{observed}, or \code{none}.
#'    Defaults to \code{expected}. However, if expected Fisher Information cannot be
#'    computed, the Monte Carlo method is used.
#' @param montecarlo_samples The number of simulations to be used for Monte Carlo
#'   Fisher Information. Defaults to 1e5.
#' @param verbose Logical. Indicating whether an announcement is printed on the
#'   console when fitting of a new model starts.
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
){

  # Check if running on Windows
  if(!Sys.info()[["sysname"]]=="Windows") {
    stop("Sorry, but model estimation only works on Windows machines. We know this must be disappointing.")
  }



  # Input validation ----
  if (!(crit %in% c("AIC","BIC"))) {
    stop("Parameter crit must be either \"AIC\" or \"BIC\".")
  }



  # Create legacy eqn file
  eqn_file <- paste0("HMMTreeR-tmpfile-", basename(model_file))
  eqn_conv_out <- HMMTreeR:::simplify_eqn(model_file = model_file, eqn_file = eqn_file)


  # Create HMMTree data file (tab-separated, first column is dataset name)
  dataset_name <- gsub(basename(data_file), pattern = ".dat|.csv|.txt| ", replacement = "")

  dat_file <- paste0("HMMTreeR-tmpfile-", dataset_name, ".dat")

  if(grepl(x = data_file, pattern = ".csv")) {
    tmp_dat <- read.csv(file = data_file)
    if(any(grepl(x = tmp_dat[, 1], pattern = ";"))) {
      tmp_dat <- read.csv2(file = data_file)
      # print(tmp_dat)
    }
  } else {
    if(grepl(x = data_file, pattern = ".dat")) {
      tmp_dat <- read.delim(file = data_file)
    }
  }
  # print(tmp_dat)

  try_conv <- suppressWarnings(as.integer(eqn_conv_out$cat_order))
  if(any(is.na(try_conv))) {
    tmp_dat <- tmp_dat[, eqn_conv_out$cat_order]
    tmp_dat <- cbind(
      data.frame(
        id = dataset_name
      )
      , tmp_dat
    )
  }

  write.table(
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
        model = eqn_file
        , data = dat_file
        , nsubj = nsubj
        , nclass = classes
        , nruns = runs
        , fi = switch(fisher_information, "expected" = 3, "montecarlo" = 2, "observed" = 1, "none" = 0)
        , mc = montecarlo_samples
        , keep_files = FALSE
      )
    )
  } else {
    # Estimate more complex models until criterion is reached
    n_classes <- 1
    repeat {
      if(verbose) {
        if(n_classes==1) {
          cat("Fitting a latent-class MPT with one class (as a baseline).\n")
        } else {
          cat("Fitting a latent-class MPT with", n_classes, "classes.\n")
        }
      }

      tmp <- model_object(
        hmmtreec(
          model = eqn_file
          , data = dat_file
          , nsubj = nsubj
          , nclass = n_classes
          , nruns = runs
          , fi = switch(fisher_information, "expected" = 3, "montecarlo" = 2, "observed" = 1, "none" = 0)
          , mc = montecarlo_samples
          , keep_files = FALSE
        )
      )

      # Stop estimating more complex models if failcodes occur
      if(tmp$description$failcode!=0) {
        break
      }

      # Stop if negative degrees of freedom occur (seems to be a byproduct of ill-conditioned FI)
      if(any(tmp$fit_statistics[c("df_M1", "df_M2", "df_S1", "df_S2")]<0)) {
        break
      }

      # Stop estimating more complex models if criterion (i.e., AIC or BIC) increased
      if(n_classes > 1 && tmp$fit_statistics[[crit]] > res[[n_classes-1]]$fit_statistics[[crit]]) {
        break
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
  file.remove(eqn_file, dat_file)

  # return
  res
}
