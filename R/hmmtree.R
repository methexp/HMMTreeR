#' @keywords internal

hmmtreec <- function(model, data, nsubj, nclass=1, nruns=1, fi=3, mc=10000, path="./", keep_files = FALSE){

  path_to_exe <-paste0(find.package("HMMTree")[1], "/HMMTreeC.exe")
  path_ <- getwd()

  # prepare
  mc_ <- format(mc, scientific = FALSE)

  ## check if input files exist
  model_name <- gsub(model, pattern = ".eqn|.EQN", replacement = "")
  dname <- strsplit(data, "[.]")[[1]][1]
  model_ <- file.path(path_, model)
  data_ <- file.path(path_, data)
  # todo: check model and data files

  ## prepare output
  comma = 0

  out <- NULL
  outfile <- file.path(paste0(model_name, ".out"))

  if(all(file.exists(c(model_, data_)))){

    # put together parameter string & call
    pars <- paste(c(model_, data_, nsubj, nclass, nruns, fi, mc_, comma), collapse = "\n")
    control_file <- write(x = pars, file = paste0(path, "/control_file.txt"))
    system(command = paste0(path_to_exe, " ", path, "/control_file.txt"))


    # check if the call was successful and returned an output file
    # return results from file "modelfilename".out
    if(file.exists(outfile)) {
      out <- read.table(file=outfile, header=TRUE, quote="", comment.char="", row.names=NULL, stringsAsFactors=FALSE)
    }

    if(!keep_files) {
      to_remove <- intersect(
        list.files(path)
        , c(
          paste0(model_name, c(".sps", ".log", ".out"))
          , "control_file.txt"
          , "lik.out"
          , "lik.err"
          , "transfer.out"
          , "transfer.err"
          , "onetrans.out"
          , "onetrans.err"
          , "onefisch.out"
          , "onefisch.err"
          , "efischer.out"
          , "efischer.err"
          , "ofischer.out"
          , "ofischer.err"
          , "tests.out"
          , "tests.err"
        )
      )
      # print(to_remove)
      file.remove(paste0(path, "/", to_remove))
    }
  }

  return(out)
}


#' @keywords internal

model_object <- function(x) {

  object <- list()

  cols <- colnames(x)



  # model description ----

  description <- c(
    Dataset = "data"
    , Ncases = "n_subjects"
    , NParam = "n_parameters"
    , Nparam = "n_parameters"
    , FItype = "fischer_information"
    , Failcode = "failcode"
    , NoClasses = "n_classes"
  )

  object$description <- x[, intersect(cols, names(description))]
  colnames(object$description) <- description[colnames(object$description)]


  # model fit ----

  fitstats <- c(
    Lik = "log-likelihood"
    , AIC = "AIC"
    , BIC = "BIC"
    , BIC2 = "BIC2"
    , G2 = "G-squared"
    , M1 = "M1"
    , df_M1 = "df_M1"
    , M2 = "M2"
    , df_M2 = "df_M2"
    , S1 = "S1"
    , df_S1 = "df_S1"
    , S2 = "S2"
    , df_S2 = "df_S2"
  )

  object$fit_statistics <- as.data.frame(
    matrix(NA, nrow = 1, ncol = length(fitstats))
  )
  colnames(object$fit_statistics) <- fitstats
  object$fit_statistics[, fitstats[intersect(cols, names(fitstats))]] <- x[, intersect(cols, names(fitstats))]

  # class weights ----

  weights <- cols[grepl(cols, pattern = "ClassWeight")]

  point_estimates <- weights[!grepl(weights, pattern = "lower|upper")]
  lower <- weights[grepl(weights, pattern = "lower")]
  upper <- weights[grepl(weights, pattern = "upper")]


  object$class_weights <- data.frame(
    class = as.numeric(gsub(x = point_estimates, pattern = "ClassWeight_", replacement = ""))
    , estimate = as.numeric(x[, point_estimates])
    , lower = as.numeric(x[, lower])
    , upper = as.numeric(x[, upper])
  )

  if(length(weights)==0) {
    object$class_weights <- data.frame(
      class = 1
      , estimate = 1
      , lower = NA_real_
      , upper = NA_real_
    )
  }



  # parameter estimates ----
  est_cols <- setdiff(cols, c(names(fitstats), weights, names(description)))

  test1 <- table(gsub(est_cols, pattern = "_upper|_lower", replacement = ""))
  parameter_names <- names(test1)
  par_list <- strsplit(parameter_names, "_")
  class_list <- lapply(X = par_list, FUN = function(x){x[length(x)]})
  i_class <- unlist(class_list)

  for (i in 1:length(parameter_names)) {
    parameter_names[i] <- gsub(parameter_names[i], pattern = paste0("_", i_class[i]), replacement = "")
  }
  # print(parameter_names)

  object$parameter_estimates <- data.frame(
    parameter = parameter_names
    , class = i_class
    , estimate = as.numeric(x[, paste0(parameter_names, "_", i_class)])
    , lower = as.numeric(x[, paste0(parameter_names, "_", i_class, "_lower")])
    , upper = as.numeric(x[, paste0(parameter_names, "_", i_class, "_upper")])
  )
  # print(paste0(parameter_names, "_", i_class))

  # return
  class(object) <- c("lc_mpt", class(object))
  object
}


#' Fit latent-class MPT models
#'
#' This function estimates one or more latent-class MPT models.
#' Starting from a one-class solution, it increments the number of latent classes
#' until a criterion is met.
#'
#' @references
#'   Stahl, C., & Klauer, K.C. (2007). HMMTree: A computer program for hierarchical multinomial processing tree models. \emph{Behavior Research Methods}, \emph{39}, 267-273.
#'
#' @export

lc <- function(model, data, nsubj, nclass_max=5, nruns=1, fi=3, mc=10000, comma=0, path="./", crit="AIC"){


  # Input validation ----
  if (!(crit %in% c("AIC","BIC"))) {
    stop("Parameter crit must be either \"AIC\" or \"BIC\".")
  }

  res <- list()

  # test fit: is (co)variance accounted for by single-class model?
  #pval <- pchisq(q=as.numeric(res[[nclass]]$S1), df=as.numeric(res[[nclass]]$S1_df), lower.tail=FALSE)
  #if ((length(pval)>0) || (pval < .05)){ # aggregate does not fit, estimate latent-class models

  n_classes <- 1
  repeat {
    res[[n_classes]] <- model_object(
      hmmtreec(
        model = model
        , data = data
        , nsubj = nsubj
        , nclass = n_classes
        , nruns = nruns
        , fi = fi
        , mc = mc
        , path=path
        , keep_files = FALSE
      )
    )

    # Stop estimating more complex models if criterion (i.e., AIC or BIC) increased
    if(n_classes > 1 && res[[n_classes]]$fit_statistics[[crit]] > res[[n_classes-1]]$fit_statistics[[crit]]) {
      break
    }
    n_classes <- n_classes + 1
    if(n_classes > nclass_max) {
      break
    }
  }

  # return
  class(res) <- c("lc_mpt_list", class(res))
  res
}


#' Fit Statistics for Latent-class MPT models
#'
#' Extracts fit statistics from latent-class MPT model objects.
#'
#' @export

fit_statistics <- function(x, ...) {
  UseMethod("fit_statistics")
}

#' @export

fit_statistics.lc_mpt <- function(x, ...) {
  x$fit_statistics
}

#' @importFrom dplyr bind_rows
#' @export

fit_statistics.lc_mpt_list <- function(x, ...) {
  dplyr::bind_rows(lapply(X = x, FUN = fit_statistics))
}


#' Parameter Estimates from Latent-class MPT models
#'
#' Extracts parameter estimates from latent-class MPT model objects.
#'
#' @rdname parameter_estimates
#' @export

parameter_estimates <- function(x, ...) {
  UseMethod("parameter_estimates")
}

#' @rdname parameter_estimates
#' @export

parameter_estimates.lc_mpt <- function(x, ...) {
  x$parameter_estimates
}

#' @rdname parameter_estimates
#' @export

parameter_estimates.lc_mpt_list <- function(x, ...) {
  lapply(x, parameter_estimates)
}


#' @rdname weighted_means
#' @export

weighted_means <- function(x, ...) {
  UseMethod("weighted_means")
}

#' @rdname weighted_means
#' @export

weighted_means.lc_mpt <- function(x, ...) {

  # Calculate variances of parameter estimates
  # CS: variances = (CIs/1.96)^2
  x$parameter_estimates$variance <- ((x$parameter_estimates$upper - x$parameter_estimates$estimate)/qnorm(p = 0.975))^2

  class_weights <- x$class_weights$estimate
  names(class_weights) <- x$class_weights$class

  x$parameter_estimates$weights <- class_weights[as.character(x$parameter_estimates$class)]
  x$parameter_estimates$wm <- x$parameter_estimates$estimate *x$parameter_estimates$weights
  x$parameter_estimates$wv <- x$parameter_estimates$variance *x$parameter_estimates$weights

  agg <- aggregate(formula = cbind(wm, wv) ~ parameter, data = x$parameter_estimates, FUN = sum)
  ci <- sqrt(agg$wv) * qnorm(p = 0.975)

  out <- data.frame(
    parameter = agg$parameter
    , estimate = agg$wm
    , lower = agg$wm - ci
    , upper = agg$wm + ci
  )

  # return
  out
}

#' @rdname weighted_means
#' @export
#'
weighted_means.lc_mpt_list <- function(x, ...) {
  lapply(X = x, weighted_means)
}




# lc_pars_wm <- function(model_list){
#
#   # select second.to.last model
#   model <- model_list[[length(model_list)-1]]
#   # how many classes & params has that model?
#   noclass <- model$NoClasses
#   noparam <- model$NParam - noclass + 1
#   nocore <- noparam/noclass
#   # class weights (row1: mean, row2: lower, row3: upper)
#   cw <- matrix(as.numeric(model[1,19:(18+noclass*3)]), nrow=3)
#   # parameter estimates (row1: mean, row2: lower, row3: upper)
#   pe <- matrix(as.numeric(model[1,(19+noclass*3):(18+noclass*3+noparam*3)]), nrow=3)
#   # only the means (cols: core model pars; rows: latent classes)
#   pem <- matrix(pe[1,], nrow=noclass)
#   # variances = (CIs/1.96)^2
#   var <- matrix( ((pe[3,]-pe[2,])/1.96)^2, nrow=noclass)
#
#   # compute weighted means
#   wm <- 1:nocore
#   for (i in 1:nocore){ # for each core model par
#     wm[i] <- sum(cw[1,]*pem[,i])
#   }
#   # compute CIs from weighted/pooled SEs
#   wci <- 1:nocore
#   for (i in 1:nocore){ # for each core model par
#     wci[i] <- 1.96*sqrt(sum(cw[1,]*var[,1])) # SE is sqrt of weighted sum of variances
#   }
#   # compute CIs from CI lengths
#   out <- rbind(weighted_mean=wm, CI_lower=wm-wci/2, CI_upper=wm+wci/2)
#   # par names
#   pnames <- names(model)[(19+noclass*3):(18+noclass*3+noparam*3)]
#   pn <- matrix(pnames, nrow=3)
#   pc <- matrix(pn[1,], nrow=3)
#   dn <- 1:nocore
#   for(i in 1:nocore){
#     dn[i] <- strsplit(pc[1,], "_")[[i]][1]
#   }
#   dimnames(out)[[2]] <- dn
#
#  return(out)
# }



