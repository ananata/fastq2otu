#' Learn Sequence Errors
#' Enables exection of DADA2's learnErrors functions. 
#' Multithread and verbose parameters are true by default
#' @param file path to FASTQ file
#' @return error object
#' @export
learnSeqErrors <- function(file) {
  # Execute the function
  errF <- dada2::learnErrors(file, multithread = TRUE, verbose = TRUE)

  # Return derep object
  return(errF)
}

