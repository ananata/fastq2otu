makeSeqsTable <- function(dadaObj, sample.name) {
  # Load required libraries
  # Extract required parameters
  bandSize <- options$dadaBandSize
  omgA <- options$dadaOmegaA

  # Execute the function
  seq.tab <- dada2::makeSequenceTable(dadaObj, orderBy = "abundance")
  rownames(seq.tab) <- sample.name

  # Return derep object
  return(seq.tab)

}

