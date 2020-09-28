#' Remove primers from sequences 
#'
#' Uses DADA2's removePrimers to remove adapter sequences from FASTQ reads.
#' To use the function only a single parameter is required. However, to use the bbduk script correctly the 
#' fastq2otu object must have its listOfAdapters, listOfSamples, and pathToRawFastq attributes defined. 
#' The listOfAdapters should define the path to a text file containing all adapter sequences to remove. 
#' The listOfSamples should define the path to a text file containing all SRA sample IDs 
#' By default, the trimmed/output files will be written to the path defined in the inDir attribute of the fastq2otu class object.
#' The path selected for the output files may or may not exist prior to the execution of this function.
#' To learn more about the fastq2out class and its attributes, please refer back to the documentation.
#' Creates output directory if it does not already exist 
#' @param object An S4 object of fastPrimerTrim
#' @return file path(s) to trimmed data
#' @export
trimAdapters <- function(object) {
	adapters <- base::readLines(con = object@listOfAdapters)
	if (object@isPaired) {
	  if (length(adapters) == 2) {
	    fAdapter <- adapter[1]
	    rAdapter <- adapter[2]
	    dada2::removePrimers(fn = object@pathToRawFastq, fout = object@pathToNoPrimers, 
			primer.fwd = fAdapter, primer.rev = rAdapter, max.mismatch = object@maxMismatch,
			allow.indels = object@allowIndels, trim.fwd = object@trimForwardReads, trim.rev = object@trimReverseReads,
			orient = object@orientReads, compress = object@compressFiles, verbose = object@verbose)
	    return(object@pathToNoPrimers)
	  } else {
	    stop(paste0("Expected two adapter sequences from file, got ", length(adapters), " instead."))
	  }
	} else {
	  if (length(adapters) == 1) {
	    dada2::removePrimers(fn = object@pathToRawFastq, fout = object@pathToNoPrimers, 
		primer.fwd = adapters, max.mismatch = object@maxMismatch, allow.indels = object@allowIndels, 
		trim.fwd = object@trimForwardReads, trim.rev = object@trimReverseReads,
                orient = object@orientReads, compress = object@compressFiles, verbose = object@verbose)
	    return(object@pathToNoPrimers)
	  } else {
	     	stop(paste0("Expected one adapter sequence from file, got ", length(adapters), " instead."))
	  }
	}		
}

