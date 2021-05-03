#' Filter and trim
#' Filter and trim paired or single-end data
#' @param object S4 object of type class fastq2otu (i.e. fastq2otu_single, fastq2otu_paired)
#' @param sample.names character vector containing sample ids
#' @param forwardFs character vecter containing paths to forward/single reads
#' @param reverseRs character vector containing paths to reverse reads
#" @importFrom dada2 filterAndTrim
#' @export
filtTrim <- function(object, sample.names, forwardFs = NA, reverseRs = NA) {
	# Extract required inputs
	label <- object@projectPrefix
	path <- object@outDir
	
	if (object@isPaired) {
		# Create output file names for filtered sequences
		filtFs <- file.path(path, paste0(label, "_filtered"), paste0(sample.names, "_R1_filt_trimmed.fastq.gz"))
		undownloaded_filtFs <- filtFs[lapply(filtFs, file.exists) == FALSE]
		forwardFs <- forwardFs[undownloaded_filtFs]
		filtRs <- file.path(path, paste0(label, "_filtered"), paste0(sample.names, "_R2_filt_trimmed.fastq.gz"))
		undownloaded_filtRs <- filtRs[lapply(filtRs, file.exists) == FALSE]
		reverseRs <- reverseRs[undownloaded_filtRs]
		if (length(undownloaded_filtFs) != 0) { 
			# Filter and Trim (only executes paths that do not already exist)
			filt_trim <- dada2::filterAndTrim(fwd=forwardFs, filt=undownloaded_filtFs, rev=reverseRs, filt.rev=undownloaded_filtRs, 
					maxEE=as.vector(object@filtMaxEE), trimLeft=as.vector(object@filtTrimLeft), trimRight=as.vector(object@filtTrimRight), 
					truncLen=as.vector(object@filtTruncLen), multithread=object@filtMultiThread,verbose=object@filtVerbose, 
					minLen=as.vector(object@filtMinLen), matchIDs=object@filtMatchIDs, truncQ=as.vector(object@filtTruncQ),
					compress=TRUE)
		
			# Save important R objects (needed for summary table)
			saveFilt <- filt_trim
			save(saveFilt, file = file.path(path, "filtered_objects.RData"))
			
			# Return both paths
			return(c(filtFs, filtRs))
		}			

	} else {		
		# Create output file names for filtered sequences
		filtFs <- file.path(path, paste0(label, "_filtered"), paste0(sample.names, "_filt_trimmed.fastq"))
		undownloaded_filtFs <- filtFs[lapply(filtFs, file.exists) == FALSE]
		forwardFs <- forwardFs[undownloaded_filtFs]
		print(length(undownloaded_filtFs) == length(forwardFs))

		if (length(undownloaded_filtFs) != 0) {
			# Filter and Trim (only executes filtered paths that do not already exist)
			filtTrim <- dada2::filterAndTrim(fwd=forwardFs, filt=undownloaded_filtFs, maxEE=object@filtMaxEE, trimLeft=object@filtTrimLeft,
								trimRight=object@filtTrimRight, truncLen=object@filtTruncLen, multithread=object@filtMultiThread,
								verbose=object@filtVerbose, minLen=object@filtMinLen, matchIDs=object@filtMatchIDs, truncQ=object@filtTruncQ)
		
			# Save important R objects (needed for summary table)
			saveFilt <- filtTrim
			save(saveFilt, file = file.path(path, "filtered_objects.RData"))
		}
	}

	# Return paths
	return(filtFs)
}

