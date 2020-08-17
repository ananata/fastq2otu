#' Import all methods
#' @import methods
NULL


#' Execute project workflow from FASTQ input to OTU output using DADA2 Workflow
#' 
#' The analysis takes place in multiple steps beginning at the creation of a central output directory. Once the directory
#' is created, a log-file is initilized that will contain all messages produced by DADA2. Accompanying the log file is a summary 
#' table documenting changes in read frequency following filtering and trimming procedures. This summary table can be used 
#' as a reference when determining the best way to modify different parameters input into the different functions.
#'
#' @param plotQuality Default is TRUE. Sequence quality distribution is plotted using DADA2's plotQualityProfile method.
#' @param isPaired Default if FALSE. If TRUE workflow for paired end data is executed. 
#' @param mergeSamples Default is TRUE. If FALSE, sample OTU tables will not be merged across.
#' @param configFile Path to config file (YML-formatted). 
#' @param downloadSeqs Default is FALSE. If TRUE, users will be able to retrieve FASTQ files from SRA using the fastq2otu's getSeqs method. 
#' @param trimAdapters Default is FALSE. If TRUE, users will be able to remove adapters sequences from data using BBTools bbduk.sh script.
#' @param generateReport Default is FALSE. If TRUE, allows a FASTQC report to be generated from input data.
#' 
#' @export
runPipeline <- function(configFile, isPaired = FALSE, plotQuality = TRUE, mergeSamples = TRUE, downloadSeqs = FALSE,
				trimAdapters = FALSE, generateReport = FALSE) {
	if (!file.exists(configFile)) {
			stop(sprintf("'%s' could not be found, please enter a valid path", configFile))
	}
	# Load all methods in package
	devtools::load_all()
	
	# Parse YAML file
	options <- yaml::yaml.load_file(configFile)
		
	# Get path to input directory (if valid)
	fp <- options$pathToData

	# Get and set path to output directory
	out <- options$outDir
	
	# Get current working directory
	curr.wd <- getwd()
	
	# Change working directory to output directory
	path <- out
	base::setwd(path)
	
	# Save all outputs to file
	if (!is.null(options$projectPrefix)) {
		con <- file(paste0(options$projectPrefix, "_fastq2otu_output.log"))
	} else {
		con <- file("fastq2otu_output.log")
		options$projectPrefix <- "myproject" # Changed to default
	}
	
	sink(con, append=TRUE)
	sink(con, append=TRUE, type="message")

		# Print the date
		message("Date: ", Sys.Date())
		message("\nAnalyzing data for ", options$projectPrefix)

		# Print package versions
		message("R version: ", version$version.string)
		message("DADA2 Version: ", packageVersion("dada2"))
		message("YAML Version: ", packageVersion("yaml"))
		if (generateReport) { message("FASTQCR Version: ", packageVersion("fastqcr")) }
		
		# Download sequences
		if (downloadSeqs) {
			object <- readConfig(configFile, type = "seqdump")
			fp <- getSeqs(object)
		}
		
		# Remove primers and update path
		if (trimAdapters) {
			message("==== Removing Primers ====")
			object <- readConfig(configFile, type = "primertrim")
			fp <- removePrimers(object)
		} else {
			fp <- options$inDir
		}
		
		# Run FASTQCR
		if (generateReport & !file.exists(file.path(options$pathToFastqcResults, paste0(options$projectPrefix, "_fastqc_report.html")))) {
			message("==== Generating FASTQC Files for Data====")
			object <- readConfig(configFile, type = "report")
			results <- runFastqcr(object, fp)
		}

		# Get file extension pattern
		if (!is.null(options$fastaPattern)) {
			REGEX_PAT <- options$fastaPattern
		} else {
			REGEX_PAT <- "^.*[1,2]?.fastq(.gz)?$"
		}
			
		if (isPaired) {
			# Extract SRA sample ids
			sample.names <- sapply(strsplit(basename(sort(list.files(fp, pattern = REGEX_PAT, full.names = TRUE))), REGEX_PAT), `[`, 1)
			
			# === Analyze as paired-end data ===
			amplicons <- paired_analysis(fp, sample.names, configFile)
			
			# Generate read frequency tables
			freq.seqtabs <- mapply(makeSeqsTable, dadaObj=amplicons, sample.name = sample.names, SIMPLIFY = FALSE)
			
		} else {
			# Extract SRA sample ids
			sample.names <- sapply(strsplit(basename(sort(list.files(fp, pattern = REGEX_PAT, full.names = TRUE))), REGEX_PAT), `[`, 1)
			
			# === Analyze as single-end data ===
			amplicons <- single_analysis(fp, sample.names, configFile, REGEX_PAT)
			
			# Generate read frequency tables
			freq.seqtabs <- mapply(makeSeqsTable, dadaObj=amplicons, sample.name=sample.names, SIMPLIFY = FALSE)
		}
		
		# Remove or Label Chimeric Sequences
		message("==== Finding Chimeric Sequences ====")
		# Get sequences
		seqtabs.nochim <- mapply(removeChimeras, seqtab=freq.seqtabs, object=object, SIMPLIFY = FALSE)

		# Save Sequence Tables
		message("==== Saving Sequence Tables ====")
		message("If the option to create a Chimera Detection Table is true (refer to config file), then two versions of the sequence table will be saved.")
		if (as.logical(options$createChimeraDetectionTable) == TRUE) {
			saveSeqs(seqtabs.nochim, sample.names, options$outDir) # Slightly modified file name will be used
			seqtabs.nochim <- freq.seqtabs # Replace with integer table
		}
		pathToSeqTables <- saveSeqs(seqtabs.nochim, sample.names)

		# Assign Taxonomy
		message("==== Assigning Taxa ====")
		otuTabs <- mapply(assignSeqTaxonomy, seqtab=freq.seqtabs, object=object, SIMPLIFY = FALSE)

		# Save OTU Tables
		message("==== Saving OTU Tables ====")
		pathToOTUTables <- saveTaxonomyTables(otuTabs, sample.names)
		
		# Merge tables
		if (options$mergeSamples) {
			message("==== Merging Tables ====")
			finalTable <- mergeSamples(pathToOTUTables, pathToSeqTables, options$finalMergedTable)
		}

		if (typeof(finalTable) == "character") {
		  message("Final table was not created")
		}

		# Restore output to console
		sink()
		sink(type="message")
		
		# Switch working directory back to original
		setwd(curr.wd)
	
		print("Done!")
}

single_analysis <- function(fp, sample.names, file, pattern = "^.*[1,2]?.fastq(.gz)?$") {
	# Sort path to extract all .fastq files (external variables can be accessed)
	Fs <- sort(list.files(fp, pattern = pattern, full.names = TRUE))

	# Filter and Trim (generates "filtered_objects.RData" file)
	message("==== Filtering and Trimming Amplicon Sequences ====")
	filtFs <- filtTrim(Fs, object) #TODO: Allow object to be passed to filtTrim function 

	# Create object - simplifies debugging process
	object <- readConfig(file, isPaired = FALSE, type = c('auto', 'filter'))
	
	# Dereplicate Sequences
	message("==== Dereplicating sequences ====")
	derepFs <- lapply(filtFs, dada2::derepFastq, derepN = object@derepN)

	# Learn Errors
	message("==== Learning Errors ====")
	errFs <- lapply(filtFs, learnSeqErrors)

	# Plot Errors
	if (!is.null(object@saveErrorsPlotPDF) & object@saveErrorsPlotPDF) {
		message("==== Plotting Errors ====")
		pdf("learned_errors_plot.pdf")
			lapply(errFs, dada2::plotErrors)
		dev.off()
	}
	
	# Denoise Data
	message("==== Removing Learned Errors ====")
	dadaFs <- mapply(dadaSeqs, derep=derepFs, err=errFs, object=object, SIMPLIFY = FALSE)
	
	# Track changes
	getSeqN <- function(x) sum(getUniques(x))
	if (nzchar(system.file(package = "HTSeqGenie"))) {
		read.in <- lapply(Fs, HTSeqGenie::getNumberOfReadsInFASTQFile)
		read.out <- lapply(filtFs, HTSeqGenie::getNumberOfReadsInFASTQFile)
		track <- cbind(sample.names, read.in, read.out, sapply(derepFs, getSeqN), sapply(dadaFs, getSeqN))
	} else {
		load("filtered_objects.RData") # Loads saveFilt object to current environment
		track <- cbind(saveFilt, sapply(derepFs, getSeqN), sapply(dadaFs, getSeqN))
	}
	
	# Specify column names
	colnames(track)[4] <- "derep"
	colnames(track)[5] <- "denoised"

	# Show tracking and save to file
	s.print <- paste0(options$projectPrefix, "_read_retention_table.txt")
	write.table(track, file = s.print, sep = "\t")

	return(dadaFs)
}

paired_analysis <- function(fp, sample.names, object) {
	# Sort path to extract all .fastq files (allows .fastq, .fq or .fastq.gz file extensions)
	Fs <- sort(list.files(fp, pattern = "*_1.fastq(.gz)?$", full.names = TRUE))
	Rs <- sort(list.files(fp, pattern = "*_2.fastq(.gz)?$", full.names = TRUE))
	
	# Filter and Trim
	message("==== Filtering and Trimming Paired-End Amplicon Sequences ====")
	filtered.files <- mapply(filtTrim, object=object, sample.names=sample.names)

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
	errFs <- lapply(filt.forward, dada2::learnErrors)
	errRs <- lapply(filt.reverse, dada2::learnErrors)

	# Plot Errors
	message("==== Plotting Errors ====")
	pdf(options$errPDF)
	  lapply(errFs, plotSeqErrors)
	  lapply(errRs, plotSeqErrors)
	dev.off()

	# Denoise Data
	message("==== Removing Learned Errors ====")
	dadaFs <- mapply(dadaSeqs, derep=derepFs, err=errFs, object=object, SIMPLIFY = FALSE)
	dadaRs <- mapply(dadaSeqs, derep=derepRs, err=errRs, object=object, SIMPLIFY = FALSE)

	# Merge forward and reverse reads
	merged_amplicons <- mapply(mergeSeqPairs, dadaFS=dadaFs, dadaRS=dadaRs, derepFS=derepFs, derepRS=derepRs, object=object, SIMPLIFY = FALSE)
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
	
	forward.derep <- sapply(derepFs, getSeqN)
	reverse.derep <- sapply(derepRs, getSeqN)
	forward.dada <- sapply(dadaFs, getSeqN)
	reverse.dada <- sapply(dadaRs, getSeqN)
	merged <- sapply(merged_amplicons, getSeqN)

	if (nzchar(system.file(package = "HTSeqGenie"))) {
		forward.in <- lapply(Fs, HTSeqGenie::getNumberOfReadsInFASTQFile)
		forward.out <- lapply(filt.forward, HTSeqGenie::getNumberOfReadsInFASTQFile)
		reverse.in <- lapply(Rs, HTSeqGenie::getNumberOfReadsInFASTQFile)
		reverse.out <- lapply(filt.reverse, HTSeqGenie::getNumberOfReadsInFASTQFile)
			
		# Combine all information
		track <- cbind(sample.names, forward.in, forward.out, reverse.in, reverse.out, forward.derep, reverse.derep, forward.dada, reverse.dada, merged)
	} else {
		load("filtered_objects") # Load saveFilt obect in current environment
		
		# Combine all information
		track <- cbind(sample.names, forward.in, forward.out, reverse.in, reverse.out, forward.derep, reverse.derep, forward.dada, reverse.dada, merged)
	}

	# Show tracking and save to file
	s.print <- paste0(options$projectPrefix, "_", options$finalSummaryTable, ".txt")
	write.table(track, file = s.print, sep = "\t")

	return(merged_amplicons)	
}


