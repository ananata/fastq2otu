#' Generate Sequence Frequency Table
#' Uses DADA2 makeSequenceTable function to generate a sequence (ASV) frequency table
#' that provides the number of reads per ASV. 
#' @param dadaObj Object returned from mergePairs or dada functions
#' @return matrix with 1-2 rows. ASVs are listed as column names and each row contains frequency information belonging to a single sample.  
#' @export
makeSeqsTable <- function(dadaObj, sample.name) {
  # Execute the function
  seq.tab <- dada2::makeSequenceTable(dadaObj, orderBy = "abundance")
  rownames(seq.tab) <- sample.name

  # Return derep object
  return(seq.tab)

}

