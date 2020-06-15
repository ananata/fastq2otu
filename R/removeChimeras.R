removeChimeras <- function(seqtab) {
  # Load required libraries
  # Determine if chimera detection table is created
  # (i.e. Identified chimeric sequences are not removed but are labeled)
  # Allows user to be able to review flagged sequences to verify chimeric characteristics
  if (options$createChimeraDetectionTable == TRUE) {
    # Extract required parameters
    minSampleFraction <- options$chimeraDetectionMinSampleFraction
    ignoreNNegatives <- options$chimeraDetectionIgnoreNegatives
    minFoldParentOverAbundance <- options$chimeraDetectionMinFoldParentOverabundance
    minParentAbundance <- options$chimeraDetectionParentAbundance
    allowOneOff <- options$chimeraDetectionAllowOneOff
    minOneOffParentDistance <- options$chimeraDetectionMinOneOffParentDistance
    maxShift <- options$chimeraDetectionMaxShift
    multithread <- options$chimeraDetectionMultiThread
    verbose <- options$chimeraDetectionVerbose

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
    verbose <- options$chimeraDetectionVerbose

    # Execute function
    seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, method="consensus", verbose=TRUE)

    # Return table with chimeric sequences removed
    return(seqtab.nochim)
  }

}

