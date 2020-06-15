assignSeqTaxonomy <- function(seqtab) {
  # Load required libraries
  # Extract required parameters
  reference <- options$taxDatabase
  minBootstrap <- options$assignTaxMinBootstrap
  tryComplement <- options$assignTaxTryComplement
  showBootstraps <- options$assignTaxOutputBootstraps
  taxLevels <- options$assignTaxLevels
  multithread <- options$assignTaxMultiThread
  verbose <- options$assignTaxVerbose

  # Execute assignTaxonomy()
  otu.tab <- dada2::assignTaxonomy(as.matrix(seqtab), reference, minBoot = minBootstrap,
                       tryRC = tryComplement, outputBootstraps = showBootstraps,
                       taxLevels = taxLevels, multithread = multithread,
                       verbose = verbose)

  # Return OTU table
  return(otu.tab)

}

