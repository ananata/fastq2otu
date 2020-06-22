#' Remove Chimeric Sequences
#' Enables the use of DADA2's isBimeraDenovoTable and/or removeBimeraDenovo functions. 
#' @param seqtab Frequency table(s) returned from makeSeqsTable function
#' @export
removeChimeras <- function(seqtab, object) {
  # Load required libraries
  # Determine if chimera detection table is created
  # (i.e. Identified chimeric sequences are not removed but are labeled)
  # Allows user to be able to review flagged sequences to verify chimeric characteristics
  if (object@createChimeraDetectionTable == TRUE) {
    # Extract required parameters
    minSampleFraction <- object@chimeraDetectionMinSampleFraction
    ignoreNNegatives <- object@chimeraDetectionIgnoreNegatives
    minFoldParentOverAbundance <- object@chimeraDetectionMinFoldParentOverabundance
    minParentAbundance <- object@chimeraDetectionParentAbundance
    allowOneOff <- object@chimeraDetectionAllowOneOff
    minOneOffParentDistance <- object@chimeraDetectionMinOneOffParentDistance
    maxShift <- object@chimeraDetectionMaxShift
    multithread <- object@chimeraDetectionMultiThread
    verbose <- object@chimeraDetectionVerbose

    # Execute the function (simple)
    seqtab.chim <- dada2::isBimeraDenovoTable(seqtab, minSampleFraction= as.numeric(minSampleFraction),
                   ignoreNNegatives= as.numeric(ignoreNNegatives), minFoldParentOverAbundance = as.numeric(minFoldParentOverAbundance),
                   minParentAbundance = as.numeric(minParentAbundance), allowOneOff= as.logical(allowOneOff),
                   minOneOffParentDistance=as.numeric(minOneOffParentDistance), maxShift=as.numeric(maxShift),
                   multithread=as.logical(multithread), verbose=as.logical(verbose))

    # Return table with chimeric sequences labeled
    return(seqtab.chim)

  } else {
    # Extract reqired parameters
    verbose <- object@chimeraDetectionVerbose

    # Execute function
    seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, method="consensus", verbose=TRUE)

    # Return table with chimeric sequences removed
    return(seqtab.nochim)
  }

}

