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
  overhang <- object$mergePairsTrimOverhang
  overlap <- object$mergePairsMinOverlap
  mMismatch <- object$mergePairsMaxMismatch
  keepRejects <- object$mergePairsReturnRejects
  concat <- object$mergePairsJustConcatenate
  allowVerbose <- object$verbose

  message("Was able to enter function, and extract desired inputs")

  # Execute function
  if (length(dadaFS) != length(dadaRS)) { stop("dadaF and dadaR must be equal in length") }
  if (length(derepFS) != length(derepRS)) { stop("derepF and derepR must be equal in length") }

  merged.seqs <- dada2::mergePairs(dadaF=dadaFS, derepF=derepFS, dadaR=dadaRS, derepR=derepRS,
                       minOverlap=overlap, maxMismatch=mMismatch, trimOverhang=overhang,
                       returnRejects=keepRejects, justConcatenate=concat, verbose=allowVerbose)

  # Return output
  return(merged.seqs)
}


