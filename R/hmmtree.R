#' @keywords internal

hmmtreec <- function(
  model_file
  , data_file
  , nsubj
  , nclass = 1
  , nruns = 1
  , fi = 3
  , mc = 1e5
  , comma = 0
  , keep_files = FALSE
){

  path_to_exe <-paste0(find.package("HMMTreeR")[1], "/HMMTreeC.exe")

  ## prepare output
  out <- NULL
  outfile <- gsub(data_file, pattern = ".dat", replacement = ".out")

  if(!all(file.exists(c(model_file, data_file)))){
    stop("Either .eqn or .dat file could not be found.")
  }

  # put together parameter string & call
  pars <- paste(
    c(
      model_file
      , data_file
      , nsubj
      , nclass
      , nruns
      , fi
      , format(mc, scientific = FALSE)
      , comma
    )
    , collapse = "\n"
  )
  write(x = pars, file = "control_file.txt")

  system(
    command = paste(
      c(path_to_exe
        , file.path(getwd(), "control_file.txt")
      )
      , collapse = " "
    )
    # , show.output.on.console = TRUE
  )


  # check if the call was successful and returned an output file
  # return results from file "modelfilename".out
  out <- try(
    read.table(file=outfile, header=TRUE, quote="", comment.char="", row.names=NULL, stringsAsFactors=FALSE)
    , silent = TRUE
  )

  if(class(out) == "try-error") {
    err_files <- list.files(pattern = ".err")
    cat("Errors in HMMTreeC.exe")
    cat(unlist(sapply(err_files, FUN = readLines)), sep = "\n")
  }

  if(!keep_files) {
    to_remove <- intersect(
      sub(list.files(recursive = TRUE, full.names = TRUE), pattern = "./", replacement = "")
      , c(
        gsub(x = data_file, pattern = ".dat", replacement = ".sps")
        , gsub(x = data_file, pattern = ".dat", replacement = ".log")
        , outfile
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
    file.remove(to_remove)
  }

  if(class(out)=="try-error") {
    stop("Output from HMMTreeC.exe not found.")
  }

  return(out)
}
