#' Fit Statistics for Latent-class MPT models
#'
#' Extracts fit statistics from latent-class MPT model objects.
#'
#' @rdname fit_statistics
#' @export

fit_statistics <- function(x, ...) {
  UseMethod("fit_statistics")
}



#' @rdname fit_statistics
#' @export

fit_statistics.lc_mpt <- function(x, ...) {
  x$fit_statistics
}



#' @rdname fit_statistics
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



#' Weighted Means of Parameter Estimates
#'
#' Extract weighted means of parameter estimates.
#'
#' @param x An object of class \code{lc_mpt} or \code{lc_mpt_list}.
#'
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
  x$parameter_estimates$wv <- x$parameter_estimates$variance * x$parameter_estimates$weights^2

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
