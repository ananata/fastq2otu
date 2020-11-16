#' Filter and trim single or paired-end data
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
		filtRs <- file.path(path, paste0(label, "_filtered"), paste0(sample.names, "_R2_filt_trimmed.fastq.gz"))

		# Check to see if files exist
		if (all(file.exists(filtFs)) & all(file.exists(filtRs))) {
			return(c(filtFs, filtRs))
		} else {
			fnotExist <- filtFs[!file.exists(filtFs)]
			rnotExist <- filtRs[!file.exists(filtRs)]
			fFs <- forwardFs[!file.exists(filtFs)]
			rRs <- reverseRs[!file.exists(filtRs)]
			
			if (length(fnotExist) != length(rnotExist)) {
				stop("Unsorted filtered files.")
			}
		}
		
		# Filter and Trim (only executes paths that do not already exist)
		filt_trim <- dada2::filterAndTrim(fwd=fFs, filt=fnotExist, rev=rRs, filt.rev=rnotExist, maxEE=as.vector(object@filtMaxEE), trimLeft=as.vector(object@filtTrimLeft),
								trimRight=as.vector(object@filtTrimRight), truncLen=as.vector(object@filtTruncLen), multithread=object@filtMultiThread,
								verbose=object@filtVerbose, minLen=as.vector(object@filtMinLen), matchIDs=object@filtMatchIDs, truncQ=as.vector(object@filtTruncQ),
								compress=TRUE)
		
		# Save important R objects (needed for summary table)
		saveFilt <- filt_trim
		save(saveFilt, file = file.path(path, "filtered_objects.RData"))
		
		# Return both paths
		return(c(filtFs, filtRs))
								

	} else {		
		# Get input
		Fs <- forwardFs

		# Create output file names for filtered sequences
		filtFs <- file.path(path, paste0(label, "_filtered"), paste0(sample.names, "_filt_trimmed.fastq"))

		# Check to see if output files already exist
		if (all(file.exists(filtFs))) {
			return(filtFs)
		} else {
			notExist <- subset(filtFs, !file.exists(filtFs))
		}

		# Filter and Trim (only executes filtered paths that do not already exist)
		filtTrim <- dada2::filterAndTrim(Fs, notExist, maxEE=object@filtMaxEE, trimLeft=object@filtTrimLeft,
								trimRight=object@filtTrimRight, truncLen=object@filtTruncLen, multithread=object@filtMultiThread,
								verbose=object@filtVerbose, minLen=object@filtMinLen, matchIDs=object@filtMatchIDs, truncQ=object@filtTruncQ)
		
		# Save important R objects (needed for summary table)
		saveFilt <- filtTrim
		save(saveFilt, file = file.path(path, "filtered_objects.RData"))
	}

	# Return paths
	return(filtFs)
}

