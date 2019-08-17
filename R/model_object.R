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
    , FItype = "fisher_information"
    , Failcode = "failcode"
    , NoClasses = "n_classes"
  )

  object$description <- data.frame(
    data = NA
    , n_subjects = NA
    , n_parameters = NA
    , fisher_information = NA
    , failcode = NA
    , n_classes = 1L
  )

  for(i in cols) {
    if(i %in% names(description)) {
      object$description[[description[i]]] <- x[[i]]
    }
  }

  # if(object$description$n_classes == 1) {
  #   object$description$fisher_information <- 3
  # }
  if(object$description$fisher_information < 0) {
    object$description$fisher_information <- 0
  }
  object$description$fisher_information <- c("none", "observed", "montecarlo", "expected")[object$description$fisher_information + 1]

  # model fit ----

  fitstats <- c(
    Lik = "-2*log-likelihood"
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

  if(length(weights)==0) {
    object$class_weights <- data.frame(
      class = 1
      , estimate = 1
      , lower = NA_real_
      , upper = NA_real_
    )
  } else {

    point_estimates <- weights[!grepl(weights, pattern = "lower|upper")]

    object$class_weights <- data.frame(
      class = as.numeric(gsub(x = point_estimates, pattern = "ClassWeight_", replacement = ""))
      , estimate = as.numeric(x[, point_estimates])
      , lower = NA_real_
      , upper = NA_real_
    )

    # only extract CI if avaliable from Fisher Information
    if(any(grepl(weights, pattern = "lower"))) {
      lower <- weights[grepl(weights, pattern = "lower")]
      upper <- weights[grepl(weights, pattern = "upper")]

      object$class_weights$lower <- as.numeric(x[, lower])
      object$class_weights$upper <- as.numeric(x[, upper])
    }
  }




  # parameter estimates ----
  est_cols <- setdiff(cols, c(names(fitstats), weights, names(description)))

  test1 <- table(gsub(est_cols, pattern = "_upper|_lower", replacement = ""))
  parameter_names <- names(test1)
  par_list <- strsplit(parameter_names, "_")
  class_list <- lapply(X = par_list, FUN = function(x){x[length(x)]})
  i_class <- unlist(class_list)

  for (i in 1:length(parameter_names)) {
    parameter_names[i] <- sub(parameter_names[i], pattern = paste0("_", i_class[i]), replacement = "")
  }
  # print(parameter_names)

  object$parameter_estimates <- data.frame(
    parameter = parameter_names
    , class = i_class
    , estimate = as.numeric(x[, paste0(parameter_names, "_", i_class)])
    , lower = NA_real_
    , upper = NA_real_
  )

  # only extract CI if avaliable from Fisher Information
  if(any(grepl(est_cols, pattern = "lower"))) {
    object$parameter_estimates$lower <- as.numeric(x[, paste0(parameter_names, "_", i_class, "_lower")])
    object$parameter_estimates$upper <- as.numeric(x[, paste0(parameter_names, "_", i_class, "_upper")])
  }

  # print(paste0(parameter_names, "_", i_class))

  # return
  class(object) <- c("lc_mpt", class(object))
  object
}
