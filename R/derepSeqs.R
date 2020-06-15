# Execute function
derepSeqs <- function(file) {
  # Load required libraries
  # Extract required parameters
  isVerbose <- options$derepVerbose
  nSample <- options$derepN

  # Execute the function
  derep <- dada2::derepFastq(file, verbose=as.numeric(isVerbose), n = as.numeric(nSample))

  # Return derep object
  return(derep)
}

