#' Remove primers from sequences 
#'
#' Uses BBTools' bbduk.sh script to remove adapter sequences from FASTQ reads.
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
removePrimers <- function(object) {
	# Extract required inputs (check for valid inputs)
	adapters <- object@listOfAdapters
	check <- ifelse(file.exists(adapters), TRUE, FALSE)
	ifelse(check, check, stop(paste0(adapters, " could not be found. Please provide a valid file path to adapters.")))
	
	samples <- object@pathToSampleIDs
	check <- ifelse(file.exists(samples), TRUE, FALSE)
	ifelse(check, check, stop(paste0(samples, " could not be found. Please provide a valid file path to text file containing SRA sample IDs.")))
	
	raw_seqs <- object@pathToRawFastq
	check <- ifelse(dir.exists(raw_seqs), TRUE, FALSE)
	ifelse(check, check, stop(paste0(raw_seqs, " could not be found. Please provide a valid directory to raw fastq files.")))
	
	trimmed_seqs <- object@pathToData
	check <- ifelse(dir.exists(trimmed_seqs), TRUE, FALSE)
	ifelse(check, check, dir.create(trimmed_seqs)) # Create output directory if it does not already exist
	
	bbdukScript <- system.files('exec', 'bbduk.sh', package = 'fastq2otu')
	
	# Set system command
	command <- paste(bbdukScript, raw_seqs, trimmed_seqs, samples, adapters, sep = " ")

	# Execute command
	system(command)

	# Return path to output
	return(trimmed_seqs)
}

