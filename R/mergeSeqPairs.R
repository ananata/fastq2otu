#" Merge forward and reverse reads
#' 
#" @param dadaFS List of denoised forward reads (output of dadaSeqs function)
#" @param dadaRS List of denoised reverse reads (output of dadaSeqs function)
#" @param derepFS List of dereplicated forward reads
#" @param derepRS List of dereplicated reverse reads
#" @param object S4 object of type fastq2otu (i.e. fastq2otu_paired, fastq2otu_single)
#' @export
mergeSeqPairs <- function(dadaFS, dadaRS, derepFS, derepRS, object) {
  # Get required inputs
  overhang <- object@mergeSeqPairsTrimOverhang
  overlap <- object@mergeSeqPairsMinOverlap
  mMismatch <- object@mergeSeqPairsMaxMismatch
  keepRejects <- object@mergeSeqPairsReturnRejects
  concat <- object@mergeSeqPairsJustConcatenate
  allowVerbose <- object@mergeSeqPairsVerbose

  message("Was able to enter function, and extract desired inputs")

  # Execute function
  merged.seqs <- dada2::mergePairs(dadaF=dadaFS, derepF=derepFS, dadaR=dadaRS, derepR=derepRS,
                       minOverlap=overlap, maxMismatch=mMismatch, trimOverhang=overhang,
                       returnRejects=keepRejects, justConcatenate=concat, verbose=allowVerbose)

  # Return output
  return(merged.seqs)
}

