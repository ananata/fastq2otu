#' Filter and trim single or paired-end data
#' @param object S4 object of type class fastq2otu (i.e. fastq2otu_single, fastq2otu_paired)
#' @export
filtTrim <- function(object, sample.names) {
	# Extract required inputs
	label <- object@projectPrefix
	path <- object@outDir
	
	if (object@isPaired) {
		# Create output file names for filtered sequences
		filtFs <- file.path(path, paste0(label, "_filtered"), paste0(sample.names, "R1_filt_trimmed.fastq.gz"))
		filtRs <- file.path(path, paste0(label, "_filtered"), paste0(sample.names, "R2_filt_trimmed.fastq.gz"))

		# Check to see if files exist
		if (all(file.exists(filtFs)) & all(file.exists(filtRs))) {
			return(c(filtFs, filtRs))
		} else {
			fnotExist <- filtFs[!file.exists(filtFs)]
			rnotExist <- filtRs[!file.exists(filtRs)]
		}

		# Filter and Trim (only executes paths that do not already exist)
		filtTrim <- dada2::filterAndTrim(fwd=forwardFs, filt=fnotExist, rev=reverseFs, filt.rev=rnotExist, maxEE=as.vector(object@filtMaxEE), trimLeft=as.vector(object@filtTrimLeft),
								trimRight=as.vector(object@filtTrimRight), truncLen=as.vector(object@filtTruncLen), multithread=object@filtMultiThread,
								verbose=object@filtVerbose, minLen=as.vector(object@filtMinLen), matchIDs=object@filtMatchIDs, truncQ=as.vector(object@filtTruncQ),
								compress=TRUE)
		
		# Save important R objects (needed for summary table)
		saveFilt <- filtTrim
		save(saveFilt, file = "filtered_objects.RData")
		
		# Return both paths
		return(c(filtFs, filtRs))
								

	} else {
		# Extract SRA run numbers from filenames (i.e. SRR829719)
		sample.names <- sapply(strsplit(basename(Fs), ".fastq"), `[`, 1)

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
		save(saveFilt, file = "filtered_objects.RData")
	}

	# Return paths
	return(filtFs)
}

