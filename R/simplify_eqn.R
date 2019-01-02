#' @keywords internal

simplify_eqn <- function(model_filename, eqn_filename) {

  read_lines <- readLines(model_filename, warn = FALSE)
  read_lines <- read_lines[-1] # simply remove first line
  read_lines <- read_lines[read_lines!=""]

  #  Handling of first row
  # Some .eqn files contain number of branches, others don't
  # n_terms <- suppressWarnings(as.integer(read_lines[[1]]))
  # if(!is.na(n_terms)) {

  #}
  n_terms <- length(read_lines)

  # Remove extraneous whitespace and tabs
  repeat{
    read_lines <- gsub(read_lines, pattern = "  |\t", replacement = " ")
    if(!any(grepl(read_lines, pattern = "  |\t"))) {
      break
    }
  }
  # read_lines <- sub(read_lines, pattern = " ", replacement = ";")
  # read_lines <- sub(read_lines, pattern = " ", replacement = ";")
  read_lines <- strsplit(read_lines, split = " ")

  model <- data.frame(
    tree = unlist(lapply(X = read_lines, FUN = function(x){x[1]}))
    , cat = unlist(lapply(X = read_lines, FUN = function(x){x[2]}))
    , term = unlist(lapply(X = read_lines, FUN = function(x){x[3]}))
    , stringsAsFactors = FALSE
  )

  # print(model)
  # check if fixed parameter values are present in the model definition
  splitted_terms <- strsplit(model$term, split = "*", fixed = TRUE)
  splitted_stripped <- lapply(X = splitted_terms, FUN = gsub, pattern = "(1-", replacement = "", fixed = TRUE)
  splitted_stripped <- unlist(lapply(X = splitted_stripped, FUN = gsub, pattern = ")", replacement = "", fixed = TRUE))

  try_conversion <- suppressWarnings(as.numeric(splitted_stripped))

  if(!all(is.na(try_conversion))) {
    stop("At least one model term seems to contain a numeric constant.")
  }

  # Convert tree identifiers to integer
  model$tree <- as.integer(factor(model$tree, levels = unique(model$tree)))
  # model$tree <- model$tree[order(model$tree)]

  # Convert categories to integer, but save ordering
  tmp_cat <- model$cat
  n_cat <- length(unique(tmp_cat))
  model$cat <- as.integer(factor(model$cat))
  cat_order <- rep(NA, n_cat)
  for (i in seq_len(n_cat)) {
    cat_order[i] <- unique(tmp_cat[model$cat==i])
  }

  # Write eqn file
  writeLines(
    text = c(n_terms, paste(model$tree, model$cat, model$term, sep = " "))
    , con = eqn_filename
    , sep = "\n"
  )

  return(
    list(
      code = 0
      , cat_order = cat_order
    )
  )
}
