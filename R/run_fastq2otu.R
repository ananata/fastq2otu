#' Import all methods
#' @import methods
NULL


#' Execute project workflow from FASTQ input to OTU output using DADA2 Workflow
#' 
#' The analysis takes place in multiple steps beginning at the creation of a central output directory. Once the directory
#' is created, a log-file is initilized that will contain all messages produces by DADA2 functions as well as a summary 
#' table documenting changes in read frequency following filtering and trimming procedures. This summary table can be used 
#' as a reference when determining the best way to modify different parameters input into the different functions.
#'
#' #' @param object An S4 object of type fastq2otu. The appropriate object can be created from a config file using
#' the readConfig function (preferred method). 
#' @export
run_fastq2otu <- function(object) {
	# Get path to input directory (if valid)
	fp <- object@pathToData

	# Get and set path to output directory
	out <- object@outDir
	if (!dir.exists(out)) {
		dir.create(out)
	}
	path <- out
	base::setwd(path)
	
	# Save all outputs to file
	con <- file(paste0(object@projectPrefix, "_dada2_output.log"))
	sink(con, append=TRUE)
	sink(con, append=TRUE, type="message")

		# Print the date
		message("Date: ", Sys.Date())
		message("\nAnalyzing data for ", object@projectPrefix)

		# Print package versions
		message("R version: ", version$version.string)
		message("DADA2 Version: ", packageVersion("dada2"))
		message("YAML Version: ", packageVersion("yaml"))
		message("FASTQCR Version: ", packageVersion("fastqcr"))
		
		# Remove primers and update path
		if (object@trimPrimers) {
			message("==== Removing Primers ====")
			fp <- removePrimers()
		} else {
			fp <- object@inDir
		}
		
		# Run FASTQCR
		if (object@runFastqc & !file.exists(file.path(object@pathToFastqcResults, paste0(object@projectPrefix, "_fastqc_report.html")))) {
			message("==== Generating FASTQC Files for Data====")
			results <- fastqcr(fp)
		}

		if (object@isPaired) {
			# Extract SRA sample ids
			sample.names <- sapply(strsplit(basename(sort(list.files(fp, pattern = "*_1.fastq(.gz)?$", full.names = TRUE))), "*_1.fastq(.gz)?$"), `[`, 1)
			
			# === Analyze as paired-end data ===
			amplicons <- paired_analysis(fp, sample.names)
			
			# Generate read frequency tables
			freq.seqtabs <- mapply(makeSeqsTable, dadaObj=amplicons, sample.name=sample.names, SIMPLIFY = FALSE)
			
		} else {
			# Extract SRA sample ids
			sample.names <- sapply(strsplit(basename(sort(list.files(fp, pattern = "*.fastq(.gz)?$", full.names = TRUE))), "*.fastq(.gz)?$"), `[`, 1)
			
			# === Analyze as single-end data ===
			amplicons <- single_analysis(fp, sample.names)
			
			# Generate read frequency tables
			freq.seqtabs <- mapply(makeSeqsTable, dadaObj=amplicons, sample.name=sample.names, SIMPLIFY = FALSE)
		}
		
		# Remove or Label Chimeric Sequences
		message("==== Finding Chimeric Sequences ====")
		# Get sequences
		seqtabs.nochim <- lapply(freq.seqtabs, removeChimeras)

		# Save Sequence Tables
		message("==== Saving Sequence Tables ====")
		message("If the option to create a Chimera Detection Table is true (refer to config file), then two versions of the sequence table will be saved.")
		if (as.logical(options$createChimeraDetectionTable) == TRUE) {
			saveSeqs(seqtabs.nochim, sample.names) # Slightly modified file name will be used
			seqtabs.nochim <- freq.seqtabs # Replace with integer table
		}
		pathToSeqTables <- saveSeqs(seqtabs.nochim, sample.names)

		# Assign Taxonomy
		message("==== Assigning Taxa ====")
		otuTabs <- lapply(freq.seqtabs, assignSeqTaxonomy)

		# Save OTU Tables
		message("==== Saving OTU Tables ====")
		pathToOTUTables <- saveTaxonomyTables(otuTabs, sample.names)
		
		# Merge tables
		message("==== Merging Tables ====")
		finalTable <- mergeSamples(pathToOTUTables, pathToSeqTables)

		if (typeof(finalTable) == "character") {
		  message("Final table was not created")
		}

		# Restore output to console
		sink()
		sink(type="message")
		
		print("Done!")
}

single_analysis(fp, sample.names) {
	# Sort path to extract all .fastq files (external variables can be accessed)
	Fs <- sort(list.files(fp, pattern = "*.fastq(.gz)?$", full.names = TRUE))

	# Run FASTQCR then comment out section of code
	if (options$runFastqc & !file.exists(file.path(options$pathToFastqcResults, paste0(options$projectPrefix, "_fastqc_report.html")))) {
	message("==== Generating FASTQC Files ====")
	results <- fastqcr(fp)
	}

	# Filter and Trim (generates "filtered_objects.RData" file)
	message("==== Filtering and Trimming Amplicon Sequences ====")
	filtFs <- filtTrim(Fs)

	# Dereplicate Sequences
	message("==== Dereplicating sequences ====")
	derepFs <- lapply(filtFs, dada2::derepFastq)

	# Learn Errors
	message("==== Learning Errors ====")
	errFs <- lapply(filtFs, learnSeqErrors)

	# Plot Errors
	message("==== Plotting Errors ====")
	pdf(options$errPDF)
		lapply(errFs, plotSeqErrors)
	dev.off()
	
	# Denoise Data
	message("==== Removing Learned Errors ====")
	dadaFs <- mapply(dadaSeqs, derep=derepFs, err=errFs, SIMPLIFY = FALSE)
	
	# Track changes
	getSeqN <- function(x) sum(getUniques(x))
	read.in <- lapply(Fs, HTSeqGenie::getNumberOfReadsInFASTQFile)
	read.out <- lapply(filtFs, HTSeqGenie::getNumberOfReadsInFASTQFile)
	track <- cbind(sample.names, read.in, read.out, sapply(derepFs, getSeqN), sapply(dadaFs, getSeqN))
	
	# Specify column names
	colnames(track)[4] <- "derep"
	colnames(track)[5] <- "denoised"

	# Show tracking and save to file
	s.print <- paste0(options$projectPrefix, "_", options$finalSummaryTable, ".txt")
	write.table(track, file = s.print, sep = "\t")

	return(dadaFs)
}

paired_analysis(fp, sample.names) {
	# Sort path to extract all .fastq files (allows .fastq, .fq or .fastq.gz file extensions)
	Fs <- sort(list.files(fp, pattern = "*_1.fastq(.gz)?$", full.names = TRUE))
	Rs <- sort(list.files(fp, pattern = "*_2.fastq(.gz)?$", full.names = TRUE))
	
	# Filter and Trim
	message("==== Filtering and Trimming Paired-End Amplicon Sequences ====")
	filtered.files <- mapply(filtTrim, forwardFs=Fs, reverseFs=Rs)

	message("The filtered files")
	lapply(sort(filtered.files), message)

	# Create output file names for filtered sequences
	filt.forward <- filtered.files[grep("*_R1_filt_trimmed.fastq.gz", filtered.files)]
	filt.reverse <- filtered.files[grep("*_R2_filt_trimmed.fastq.gz", filtered.files)]

	# Dereplicate Sequences
	message("==== Dereplicating sequences ====")
	derepFs <- lapply(filt.forward, dada2::derepFastq)
	derepRs <- lapply(filt.reverse, dada2::derepFastq)

	# Learn Errors
	message("==== Learning Errors ====")
	errFs <- lapply(filt.forward, learnSeqErrors)
	errRs <- lapply(filt.reverse, learnSeqErrors)

	# Plot Errors
	message("==== Plotting Errors ====")
	pdf(options$errPDF)
	  lapply(errFs, plotSeqErrors)
	  lapply(errRs, plotSeqErrors)
	dev.off()

	# Denoise Data
	message("==== Removing Learned Errors ====")
	dadaFs <- mapply(dadaSeqs, derep=derepFs, err=errFs, SIMPLIFY = FALSE)
	dadaRs <- mapply(dadaSeqs, derep=derepRs, err=errRs, SIMPLIFY = FALSE)

	# Merge forward and reverse reads
	merged_amplicons <- mapply(mergeSeqPairs, dadaFS=dadaFs, dadaRS=dadaRs, derepFS=derepFs, derepRS=derepRs, SIMPLIFY = FALSE)
	cat(names(merged_amplicons))
	m.print <- paste0(options$projectPrefix, "_merged_pairs_table.csv")

	# A large merged table is created
	# Each table is differentiated by a prefix in the column names
	# The order in which the tables are read corresponds to the order in which they are merged
	write.table(merged_amplicons, file = m.print, sep = "\t")
	
	# Make Sequence Frequency Table
	message("==== Creating Sequence Tables ====")
	if (length(dadaFs) != length(Fs)) {
		stop("Number of dada objects does not equal number of samples")
	}
	
	# Count number of reads in sequences (load HTSeqGenie library)
	getSeqN <- function(x) sum(getUniques(x))
	forward.in <- lapply(Fs, HTSeqGenie::getNumberOfReadsInFASTQFile)
	forward.out <- lapply(filt.forward, HTSeqGenie::getNumberOfReadsInFASTQFile)
	reverse.in <- lapply(Rs, HTSeqGenie::getNumberOfReadsInFASTQFile)
	reverse.out <- lapply(filt.reverse, HTSeqGenie::getNumberOfReadsInFASTQFile)
	forward.derep <- sapply(derepFs, getSeqN)
	reverse.derep <- sapply(derepRs, getSeqN)
	forward.dada <- sapply(dadaFs, getSeqN)
	reverse.dada <- sapply(dadaRs, getSeqN)
	merged <- sapply(merged_amplicons, getSeqN)

	# Combine all information
	track <- cbind(sample.names, forward.in, forward.out, reverse.in, reverse.out, forward.derep, reverse.derep, forward.dada, reverse.dada, merged)

	# Show tracking and save to file
	s.print <- paste0(options$projectPrefix, "_", options$finalSummaryTable, ".txt")
	write.table(track, file = s.print, sep = "\t")

	return(merged_amplicons)	
}


