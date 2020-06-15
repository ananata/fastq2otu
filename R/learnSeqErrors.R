learnSeqErrors <- function(file) {
  # Load required libraries
  # Read YML file (get configFile object from main.R)
  options <- yaml.load_file(configFile)

  # Extract required parameters
  nSample <- options$errN
  nThreads <- options$errMultithread

  # Execute the function
  errF <- dada2::learnErrors(file)

  # Return derep object
  return(errF)
}

