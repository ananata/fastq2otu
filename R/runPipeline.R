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
#' @param getMergeSamples Default is TRUE. If FALSE, sample OTU tables will not be merged across.
#' @param configFile Path to config file (YML-formatted). 
#' @param downloadSeqs Default is FALSE. If TRUE, users will be able to retrieve FASTQ files from SRA using the fastq2otu's getSeqs method. 
#' @param trimAdapters Default is FALSE. If TRUE, users will be able to remove adapters sequences from data using BBTools bbduk.sh script.
#' @param generateReport Default is FALSE. If TRUE, allows a FASTQC report to be generated from input data.
#' @importFrom yaml yaml.load_file
#' @importFrom dada2 dada derepFastq learnErrors plotErrors getUniques makeSequenceTable
#' 
#' @export
runPipeline <- function(configFile, isPaired = FALSE, getQuality = TRUE, getMergedSamples = TRUE, getDownloadedSeqs = FALSE,
				getTrimmedAdapters = FALSE, getGeneratedReport = FALSE) {
	if (!file.exists(configFile)) {
			stop(sprintf("'%s' could not be found, please enter a valid path", configFile))
	}
	
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
		message("\nR version: ", version$version.string)
		message("\nDADA2 Version: ", packageVersion("dada2"))
		message("\nYAML Version: ", packageVersion("yaml"))
		if (getGeneratedReport) { message("\nFASTQCR Version: ", packageVersion("fastqcr")) }
		
		# Download sequences - tested on 8/26/2020
		if (getDownloadedSeqs == TRUE) {
			message("\n==== Downloading sequences from NCBI ====")
			object <- readConfig(configFile, type = "seqdump")
			fp <- getSeqs(object)
		} 
		
		# Remove primers and update path
		if (getTrimmedAdapters) {
			message("\n==== Removing Primers ====")
			object <- readConfig(configFile, type = "primertrim")
			fp <- trimAdapters(object)
		} else {
			fp <- options$pathToData
		}
		
		# Run FASTQCR
		if (getGeneratedReport & !file.exists(file.path(options$pathToFastqcResults, paste0(options$projectPrefix, "_fastqc_report.html")))) {

			message("\n==== Generating FASTQC Files for Data====")
			object <- readConfig(configFile, type = "report")
			results <- runFastqcr(object, fp)
		}

		# Plot Quality Distribution and save object
		if (getQuality) {
			message("\n==== Plotting quality distribution BEFORE trimming ====")
			object <- readConfig(configFile, type = "qualityplot")
			if (is.na(object)) 
				stop("Unable to generate quality plot object")

			plots <- plotQuality(fp, options$projectPrefix, object)
			
			# Write PDF files in output directory
			if (!is.null(plots)) {
				message("Created: ", paste0(options$projectPrefix, "_fastq2otu_quality_plots_BEFORE.pdf"))
				pdf(file = paste0(options$projectPrefix, "_fastq2otu_quality_plots_BEFORE.pdf"))
					plots
				dev.off()
			} else {
				stop("Unable to write quality plots to PDF") # Should rarely execute. Mostly for debugging purposes
			}
		}
			
		if (isPaired) {
			 # Get file extension pattern
                	if (!is.null(options$fastaPattern) & length(options$fastaPattern) == 2) {
                       		REGEX_PAT <- options$fastaPattern
               	 	} else {
				REGEX_PAT <- c("*_1.fastq(.gz)?", "*_2.fastq(.gz)?")
			}
			
			# === Analyze as paired-end data ===
			if (is.null(REGEX_PAT)) {
				# Extract SRA sample ids
				sample.names <- na.omit(sapply(strsplit(basename(sort(list.files(fp, pattern = "*_1.fastq(.gz)?", full.names = TRUE))),"*_1.fastq(.gz)?"), `[`, 1))
				amplicons <- paired_analysis(fp = fp, sample.names = sample.names, getQuality = getQuality, file = configFile)
			} else {
				# Extract SRA sample ids
	                        sample.names <- na.omit(sapply(strsplit(basename(sort(list.files(fp, pattern = REGEX_PAT[1], full.names = TRUE))),REGEX_PAT[1]), `[`, 1))
				amplicons <- paired_analysis(fp = fp, sample.names = sample.names, file = configFile, getQuality = getQuality, REGEX_1 = REGEX_PAT[1], REGEX_2 = REGEX_PAT[2])
			}

			# Create read frequency table
			freq.seqtabs <- dada2::makeSequenceTable(amplicons) 
			rownames(freq.seqtabs) <- sample.names
			
		} else {
			# Get file extension pattern
                        if (!is.null(options$fastaPattern) & length(options$fastaPattern) == 1) {
                                REGEX_PAT <- options$fastaPattern
                        } else {
                                REGEX_PAT <- "*.fastq(.gz)?$"
                        }
			
			# Extract SRA sample ids
			sample.names <- na.omit(sapply(strsplit(basename(sort(list.files(fp, pattern = REGEX_PAT, full.names = TRUE))), REGEX_PAT), `[`, 1))
			
			# === Analyze as single-end data ===
			amplicons <- single_analysis(fp = fp, sample.names = sample.names, file = configFile, getQuality = getQuality, REGEX = REGEX_PAT)
			
			# Generate read frequency tables
                        freq.seqtabs <- lapply(amplicons, dada2::makeSequenceTable)
			rownames(freq.seqtabs) <- sample.names

		}
		
		# Remove or Label Chimeric Sequences
		message("\n==== Finding Chimeric Sequences ====")
		# Get sequences
		seqtabs.nochim <- removeChimeras(seqtab=freq.seqtabs, object=object)		
		if (is.na(seqtabs.nochim))
			stop("Problem with seqtabs.nochim object")		

		# Save Sequence Tables
		message("\n==== Saving Sequence Tables ====")
		if (as.logical(options$createChimeraDetectionTable) == TRUE) {
			message("Output directory: ", out)

			# Creates an additional sequence table
			nothing <- saveSeqs(seqtabs.nochim, na.omit(sample.names), out, add.table = TRUE, label=options$projectPrefix) # Slightly modified file name will be used
			seqtabs.nochim <- freq.seqtabs # Replace with integer table
		}
		pathToSeqTables <- saveSeqs(seqtabs.nochim, na.omit(sample.names), out, add.table = FALSE, label=options$projectPrefix)
		message("Created: ")
		lapply(na.omit(pathToSeqTables), message)

		# Assign Taxonomy
		message("\n==== Assigning Taxa ====")
		tax.object <- readConfig(configFile, type = "assignTax")
		otuTabs <- lapply(1:length(freq.seqtabs), function(input, seqtabs, myobject) {
				return(assignSeqTaxonomy(seqtab=seqtabs, object=myobject))}, seqtabs=freq.seqtabs, myobject=tax.object)
		message("Database: ", tax.object@taxDatabase)

		# Save OTU Tables
		message("\n==== Saving OTU Tables ====")
		pathToOTUTables <- saveTaxonomyTables(otuTabs, na.omit(sample.names), out, options$projectPrefix)
		message("Created: ")
		lapply(pathToOTUTables, message)		

		# Merge tables
		if (getMergedSamples) {
			message("\n==== Merging OTU Tables ====")
			finalTable <- mergeSamples(pathToOTUTables, pathToSeqTables, options$projectPrefix)
		}

		if (typeof(finalTable) == FALSE) {
		  message("Final table was not created")
		}

		# Restore output to console
		sink()
		sink(type="message")
		
		# Switch working directory back to original working directory
		setwd(curr.wd)
	
		print("Done! All outputs were sucessfully generated.")
}

single_analysis <- function(fp, sample.names, file, getQuality = FALSE, REGEX = "*.fastq(.gz)?$") {
	# Create object - simplifies debugging process
        object <- readConfig(file, isPaired = FALSE, type = c('auto', 'filter'))

	# Sort path to extract all .fastq files (external variables can be accessed)
	Fs <- sort(list.files(fp, pattern = REGEX, full.names = TRUE))
	message("Input files: ", Fs)

	# Filter and Trim (generates "filtered_objects.RData" file)
	message("\n==== Filtering and Trimming Amplicon Sequences ====")
	if (!is.null(object)) { 
		filtFs <- filtTrim(sample.names=sample.names, object=object, forwardFs=Fs)
		message("Created: ")
		lapply(filtFs, message)
	} else {
		stop("Error created when creating fastFilt object")
	}

	# Plot Quality Distribution and save object
        if (getQuality) {
                message("\n==== Plotting quality distribution AFTER trimming ====")
                plot_object <- readConfig(file, type = "qualityplot")
                plots <- plotQuality(filtFs, object@projectPrefix, plot_object)

                # Write PDF files in output directory
                if (!is.null(plots)) {
			message("Created: ", paste0(object$projectPrefix, "_fastq2otu_quality_plots_AFTER.pdf"))
                        pdf(file = paste0(object@projectPrefix, "_fastq2otu_quality_plots_AFTER.pdf"))
                        	plots
                	dev.off()
                } else {
                	stop("Unable to write quality plots to PDF") # Should rarely execute. Mostly for debugging purposes
        	}
        }

	# Dereplicate Sequences
	message("\n==== Dereplicating sequences ====")
	derepFs <- lapply(filtFs, dada2::derepFastq)

	# Learn Errors
	message("==== Learning Errors ====")
	errFs <- lapply(filtFs, dada2::learnErrors, multithread = TRUE)

	# Plot Errors
	if (!is.null(object@saveErrorsPlotPDF) & object@saveErrorsPlotPDF) {
		message("\n==== Plotting Errors ====")
                message("Created: ", paste0(object@projectPrefix, "_learned_errors_plot.pdf"))
		pdf(paste0(object@projectPrefix, "_learned_errors_plot.pdf"))
			lapply(errFs, dada2::plotErrors, err_in=TRUE)
		dev.off()
	}
	
	# Denoise Data
        message("\n==== Removing Learned Errors ====")
        if (!is.null(object@dadaBandSize) & !is.null(object@dadaBandSize)) {
                dadaFs <- lapply(1:length(derepFs), function(input, dereps, errs) {
                                return(dada2::dada(derep=dereps[[input]], err=errs[[input]]))}, dereps = derepFs, errs=errFs)
        } else {
                dadaFs <- lapply(1:length(derepFs), function(input, dereps, errs, band_size, omega_a) {
                                return(dada2::dada(derep=dereps[[input]], err=errs[[input]]))}, dereps=derepFs, errs=errFs, band_size=object@dadaBandSize, omega_a=object@dadaOmegaA)
        }
	
	# Track changes
	getSeqN <- function(x) sum(dada2::getUniques(x))
	if (file.exists("filtered_objects.RData")) {
		load("filtered_objects.RData") # Loads saveFilt object to current environment
		track <- cbind(saveFilt, sapply(derepFs, getSeqN), sapply(dadaFs, getSeqN))
	}
	
	# Specify column names
	colnames(track)[4] <- "derep"
	colnames(track)[5] <- "denoised"

	# Show tracking and save to file
	s.print <- paste0(object@projectPrefix, "_read_retention_table.txt")
	write.table(track, file = s.print, sep = "\t")

	return(dadaFs)
}

paired_analysis <- function(fp, sample.names, file, getQuality = FALSE, REGEX_1 = "*_1.fastq(.gz)?$", REGEX_2 = "*_2.fastq(.gz)?$" ) {
	# Sort path to extract all .fastq files (allows .fastq or .fastq.gz file extensions)
	Fs <- sort(list.files(fp, pattern = REGEX_1, full.names = TRUE))
	Rs <- sort(list.files(fp, pattern = REGEX_2, full.names = TRUE))

	if (length(Fs) != length(Rs)) {
		stop("Error: Unequal number of files containing forward and reverse sequences")
	}

	# Create object - simplifies debugging process
        object <- readConfig(file, isPaired = TRUE, type = c('auto', 'filter'))

	# Filter and Trim
	message("\n==== Filtering and Trimming Paired-End Amplicon Sequences ====")
	if (!is.null(object)) {
		filtered.files <- filtTrim(sample.names=sample.names, object=object, forwardFs=Fs, reverseRs=Rs)
		message("Created: ")
		lapply(filtered.files, message)
	} else {
		stop("Error generated when creating fastFilt object")
	}

	# Create output file names for filtered sequences
	filt.forward <- filtered.files[grep("*_R1_filt_trimmed.fastq.gz", filtered.files)]
	filt.reverse <- filtered.files[grep("*_R2_filt_trimmed.fastq.gz", filtered.files)]
	
        # Plot Quality Distribution and save object
        if (getQuality) {
                message("\n==== Plotting quality distribution AFTER trimming ====")
                plot_object <- readConfig(file, type = "qualityplot")
		ordered.files <- c(rbind(filt.forward, filt.reverse)) # Alternate plots from forward and reverse sequences
                plots <- plotQuality(ordered.files, options$projectPrefix, plot_object)

                # Write PDF files in output directory
                if (!is.null(plots)) {
                        message("Created: ", paste0(object@projectPrefix, "_fastq2otu_quality_plots_AFTER.pdf"))
                        pdf(file = paste0(object@projectPrefix, "_fastq2otu_quality_plots_AFTER.pdf"))
                                plots
                        dev.off()
                } else {
                        stop("Unable to write quality plots to PDF") # Should rarely execute. Mostly for debugging purposes
                }
        }


	# Dereplicate Sequences
	message("\n==== Dereplicating sequences ====")
	derepFs <- lapply(filt.forward, dada2::derepFastq)
	derepRs <- lapply(filt.reverse, dada2::derepFastq)

	# Learn Errors
	message("\n==== Learning Errors ====")
	errFs <- lapply(filt.forward, dada2::learnErrors, multithread=TRUE)
	errRs <- lapply(filt.reverse, dada2::learnErrors, multithread=TRUE)

	# Plot Errors
	if (!is.null(object@saveErrorsPlotPDF) & object@saveErrorsPlotPDF) {
		message("\n==== Plotting Errors ====")
		message("Created: ", object@projectPrefix, "_learned_errors_plots.pdf")
		pdf(paste0(object@projectPrefix, "_learned_errors_plots.pdf"))
		  lapply(errFs, dada2::plotErrors)
		  lapply(errRs, dada2::plotErrors)
		dev.off()
	}

	# Denoise Data
	message("\n==== Removing Learned Errors ====")
	if (!is.null(object@dadaBandSize) & !is.null(object@dadaBandSize)) {
		dadaFs <- lapply(1:length(derepFs), function(input, dereps, errs) {
  				return(dada2::dada(derep=dereps[[input]], err=errs[[input]]))}, dereps = derepFs, errs=errFs)
		dadaRs <- lapply(1:length(derepRs), function(input, dereps, errs) {
                                return(dada2::dada(derep=dereps[[input]], err=errs[[input]]))}, dereps = derepRs, errs=errRs)
	} else {
		dadaFs <- lapply(1:length(derepFs), function(input, dereps, errs, band_size, omega_a) {
                                return(dada2::dada(derep=dereps[[input]], err=errs[[input]]))}, dereps=derepFs, errs=errFs, band_size=object@dadaBandSize, omega_a=object@dadaOmegaA)
                dadaRs <- lapply(1:length(derepRs), function(input, dereps, errs) {
                                return(dada2::dada(derep=dereps[[input]], err=errs[[input]]))}, dereps=derepRs, errs=errRs, band_size=object@dadaBandSize, omega_a=object@dadaOmegaA)
	}

	# Merge forward and reverse reads
        message("==== Merge Forward and Reverse Reads ====")
        if (length(dadaFs) != length(dadaRs)) {
		stop("An error occurred when denoising the data")
	}

	# Name the derep-class objects by the sample names
	names(derepFs) <- sample.names
	names(derepRs) <- sample.names

	# Merge pairs
	merged_amplicons <- mergeSeqPairs(dadaFS=dadaFs, dadaRS=dadaRs, derepFS=derepFs, derepRS=derepRs, object)
	
	message("\n==== Saving Merged Pairs Table ====") 
	m.print <- paste0(object@projectPrefix, "_merged_pairs_table.rds")
	message("Created: ", m.print)

	# A large merged table is created
	# Each table is differentiated by a prefix in the column names
	# The order in which the tables are read corresponds to the order in which they are merged
	saveRDS(merged_amplicons, file = m.print)
	
	# Make Sequence Frequency Table
	message("\n==== Creating Summary Table ====")
	if (length(dadaFs) != length(Fs)) {
		stop("Number of dada objects does not equal number of samples")
	}
	
	# Count number of reads in sequences (load HTSeqGenie library) # TODO: Find alternative
	getSeqN <- function(x) sum(dada2::getUniques(x))
	
	forward.derep <- sapply(derepFs, getSeqN)
	reverse.derep <- sapply(derepRs, getSeqN)
	forward.dada <- sapply(dadaFs, getSeqN)
	reverse.dada <- sapply(dadaRs, getSeqN)

	merged <- sapply(merged_amplicons, getSeqN)

	if (file.exists("filtered_objects.RData")) {
		load("filtered_objects.RData") # Load saveFilt obect in current environment
		
		# Combine all information
		track <- cbind(sample.names, saveFilt, forward.derep, reverse.derep, forward.dada, reverse.dada, merged)
		
	        # Show tracking and save to file
        	message("\n==== Saving Final Summary Table ====")
        	s.print <- paste0(object@projectPrefix, "_read_retention_table.txt")
        	write.table(track, file = s.print, sep = "\t")
	}

	return(merged_amplicons)	
}


