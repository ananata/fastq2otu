dadaSeqs <- function(derep, err) {
  # Load required libraries
  # Extract required parameters
  bandSize <- as.numeric(options$dadaBandSize)
  omgA <- as.numeric(options$dadaOmegaA)

  # Execute the function
  dada <- dada::dada(derep, err=err, BAND_SIZE = bandSize, OMEGA_A = omgA)

  # Return derep object
  return(dada)
}

