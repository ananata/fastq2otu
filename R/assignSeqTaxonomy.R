#' Assign Taxonomy 
#' Uses DADA2's assignTaxonomy function to assign taxa to ASV sequences. 
#'
#' @param seqtab Sequence frequency table
#' @param object S4 object of type fastq2otu (i.e. fastPaired, fastSingle)
#' @importFrom dada2 assignTaxonomy
#' @return OTU table
#' 
assignSeqTaxonomy <- function(seqtab, object) {
  # Load required libraries
  # Extract required parameters
  referenceDatabase <- object@taxDatabase
  minBootstrap <- object@assignTaxMinBootstrap
  tryComplement <- object@assignTaxTryComplement
  showBootstraps <- object@assignTaxOutputBootstraps
  taxLevels <- object@assignTaxLevels
  multithread <- object@assignTaxMultiThread
  verbose <- object@assignTaxVerbose

  # Execute assignTaxonomy()
  otu.tab <- dada2::assignTaxonomy(as.matrix(seqtab), referenceDatabase, minBoot = minBootstrap,
                       tryRC = tryComplement, outputBootstraps = showBootstraps,
                       taxLevels = taxLevels, multithread = multithread,
                       verbose = verbose)

  # Return OTU table
  return(otu.tab)

}

