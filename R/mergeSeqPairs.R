mergeSeqPairs <- function(dadaFS, dadaRS, derepFS, derepRS) {
  # Get required inputs
  overhang <- options$mergeSeqPairsTrimOverhang
  overlap <- options$mergeSeqPairsMinOverlap
  mMismatch <- options$mergeSeqPairsMaxMismatch
  keepRejects <- options$mergeSeqPairsReturnRejects
  concat <- options$mergeSeqPairsJustConcatenate
  allowVerbose <- options$mergeSeqPairsVerbose

  # Execute function
  merged <- dada2::mergePairs(dadaF = dadaFS, derepF = derepFS, dadaR = dadaRS, derepR = derepRS,
                       minOverlap = as.numeric(overlap), maxMismatch = as.numeric(mMismatch), trimOverhang = overhang,
                       returnRejects = keepRejects, justConcatenate = concat, verbose = allowVerbose)

  # Return output
  return(merged)
}

