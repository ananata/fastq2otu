# Class List
#	- fastq2otu 
# 		- fastSingle
#			- fastAssignTaxa
#			- fastFilt
# 		- fastDouble
#			- fastAssignTaxa
#			- fastFilt
#	- fastPlotQuality
#	- fastReport
#	- fastSeqDump
# 	- fastPrimerTrim
#		

#===================
# Create validation methods
#===================
#' Check for valid inputs for fastq2otu
#'
#' This function checks to see if user inputs are present and valid. Missing or invalid inputs will cause an
#' error message to be displayed.
#'
#' @param object fastq2otu-type S4 object
#' @return TRUE (boolean) or error message detailing source(s) of error.
#' @keyword internal
#' @export
check_fastq2otu <- function(object) {
	errors <- character()
		length_input <- length(object@pathToData)
	if (length_input == 0) {
		msg <- paste("No input provided for pathToData. A value is required")
		errors <- c(errors, msg)
	}
	else if (!dir.exists(object@pathToData)) {
		msg <- paste("Provided input: ", object@pathToData, " is invalid and does not exist", sep = "")
		errors <- c(errors, msg)
	}
	
    length_output <- length(object@outDir)
	if (length_output == 0) {
		msg <- paste("No input provided for outDir. A value is required")
		errors <- c(errors, msg)
	}
	else if (!dir.exists(object@outDir)) {
		msg <- paste("Provided input: ", object@outDir, " is invalid and does not exist", sep = "")
		errors <- c(errors, msg)
	}
	
	if (length(errors) == 0) TRUE else errors
}

#' Check for valid inputs for fastFilter
#'
#' This function checks to see if filtering parameters are present and valid. Missing or invalid inputs will cause an
#' error message to be displayed.
#'
#' @param object fastFilt-type S4 object
#' @return TRUE (boolean) or error message detailing source(s) of error.
#' @keyword internal
#' @export
check_filt_params <- function(object) {
	errors <- character()
	if (object@isPaired) {
		if (length(object@filtMaxEE) != 2) {
			msg <- "Invalid input. Expected two values for filtMaxEE"
			errors <- c(errors, msg)
		}
		if (length(object@filtTrucQ) != 2) { 
			msg <- "Invalid input. Expected two values for filtTruncQ"
                        errors <- c(errors, msg)
		}
		if (length(object@filtTruncLen) != 2) { 
			msg <- "Invalid input. Expected two values for filtTruncLen"
                        errors <- c(errors, msg)
		}
		if (length(object@filtTrimLeft) != 2) { 
			msg <- "Invalid input. Expected two values for filtTrimLeft"
                        errors <- c(errors, msg)
		}
		if (length(object@filtTrimRight) != 2) { 
			msg <- "Invalid input. Expected two values for filtTrimRight"
                        errors <- c(errors, msg)
		}
		if (length(object@filtMinLen) != 2) {
			msg <- "Invalid input. Expected two values for filtMinLen"
                        errors <- c(errors, msg)
		}
	}
	else {
		if (length(object@filtMaxEE) != 1) {
			msg <- "Invalid input. Expected a single value for filtMaxEE"
						errors <- c(errors, msg)
		}
		if (length(object@filtTrucQ) != 1) { 
			msg <- "Invalid input. Expected a single value for filtTruncQ"
                        errors <- c(errors, msg)
		}
		if (length(object@filtTruncLen) != 1) { 
			msg <- "Invalid input. Expected a single value for filtTruncLen"
                        errors <- c(errors, msg)
		}
		if (length(object@filtTrimLeft) != 1) { 
			msg <- "Invalid input. Expected a single value for filtTrimLeft"
                        errors <- c(errors, msg)
		}
		if (length(object@filtTrimRight) != 1) { 
			msg <- "Invalid input. Expected a single value for filtTrimRight"
                        errors <- c(errors, msg)
		}
		if (length(object@filtMinLen) != 1) {
			msg <- "Invalid input. Expected a single value for filtMinLen"
                        errors <- c(errors, msg)
		}
		
	}
	if (length(errors) == 0) TRUE else errors
}

#' Check for valid inputs for fastSeqDump
#'
#' This function checks to see if fastSeqDump parameters are present and valid. Missing or invalid inputs will cause an
#' error message to be displayed.
#'
#' @param object fastSeqDump-type S4 object
#' @return TRUE (boolean) or error message detailing source(s) of error.
#' @keyword internal
#' @export
check_seq_dump <- function(object) {
	errors <- character()
	if (!file.exists(object@pathToSampleURLs) | !file.exists(object@pathToSampleIDs)) {
			if (!file.exists(object@pathToSampleURLs)) {
				msg <- paste0(object@pathToSampleURLs, " could not be found. Please enter valid file path.")
                        errors <- c(errors, msg)
			} else {
				msg <- paste0(object@pathToSampleIDs, " could not be found. Please enter valid file path.")
                        errors <- c(errors, msg)
			}
	}
	else if (!file.exists(object$pathToFastqDump) | !dir.exists(object$pathToFastqDump)) {
				msg <- "fastq-dump script could not be found."
                        errors <- c(errors, msg)
	}
	else if (is.na(object@pathToSampleURLs) & !is.na(object@pathToSampleIDs)) {
		if (!object@installFastqDump & is.na(object$pathToFastqDump)) { 
				msg <- "Please provide a path to the fastq-dump script"
                        errors <- c(errors, msg)
		}
		if (object@installFastqDump & is.na(object$pathToFastqDump)) {
				msg <- "Please specify a directory to install fastq-dump"
                        errors <- c(errors, msg)
		}
	}
	else if (is.na(object@pathToSampleURLs) & is.na(object@pathToSampleIDs)) {
				msg <- "Unable to download files. Please verify that all required parameters are supplied."
                        errors <- c(errors, msg)
	}
	
	if (length(errors) == 0) TRUE else errors
}

#' Check for valid inputs for fastAssignTaxa
#'
#' This function checks to see if fastAssignTaxa parameters are present and valid. Missing or invalid inputs will cause an
#' error message to be displayed.
#'
#' @param object fastSeqDump-type S4 object
#' @return TRUE (boolean) or error message detailing source(s) of error.
#' @keyword internal
#' @export
check_assign_tax <- function(object) {
	errors <- c()
	if (!file.exists(object@taxDatabase)) {
		msg <- paste("Provided input: ", object@taxDatabase, " is invalid and does not exist", sep = "")
		errors <- c(errors, msg)
	}
	if (length(errors) == 0) TRUE else errors
}


# Can extend and create with new("fastq2otu", ...)
#' FASTQ2OTU Custom Super-Class
#'
#' Use to convert parameters provided by user to accessible attributes associated with
#' an S4-type data class.
#' @keyword internal
#' @export
setClass("fastq2otu",
                slots = c(
					# === Required Inputs ===
					projectPrefix = "character",
					outDir = "character",
					isPaired = "logical", 
					
					# === Provide other options ===
					downloadSeqs = "logical", 
					trimAdapters = "logical",
					generateReport = "logical", 
					mergeSamples = "logical",
					plotQuality = "logical",
					
					# === Dereplicate reads to keep only unique sequences ===
					derepVerbose = "logical",
					derepN = "numeric",

					# === Detect and learn error patterns in sequences ===
					errN = "numeric",
					errMultithread = "logical",
					saveErrorsPlotPDF = "logical",

					# === Denoise data ===
					dadaBandSize = "numeric",
					dadaOmegaA = "numeric",
					
					# === Find chimeric sequences ===
					createChimeraDetectionTable = "logical",
					chimeraDetectionMinSampleFraction = "numeric",
					chimeraDetectionIgnoreNegatives = "numeric",
					chimeraDetectionMinFoldParentOverabundance = "numeric",
					chimeraDetectionParentAbundance = "numeric",
					chimeraDetectionAllowOneOff = "logical",
					chimeraDetectionMaxShift = "numeric",
					chimeraDetectionMultiThread = "logical",
					chimeraDetectionVerbose = "logical"
			),
		prototype = list(
				projectPrefix = "myproject",
				outDir = "NA_character_",
				isPaired = FALSE, 

				downloadSeqs = FALSE, 
				trimAdapters = FALSE,
				generateReport = FALSE,
				mergeSamples = FALSE,
				plotQuality = FALSE, 
					
				derepVerbose = TRUE,
				derepN = 1e+06,

				errN = 1e+08,
				errMultithread = FALSE,
				saveErrorsPlotPDF = FALSE,

				dadaBandSize = 16,
				dadaOmegaA = 1e-40,

				createChimeraDetectionTable = FALSE,
				chimeraDetectionMinSampleFraction = 0.9,
				chimeraDetectionIgnoreNegatives = 1,
				chimeraDetectionMinFoldParentOverabundance = 1.5,
				chimeraDetectionParentAbundance = 2,
				chimeraDetectionAllowOneOff = FALSE,
				chimeraDetectionMaxShift = 16,
				chimeraDetectionMultiThread = TRUE,
				chimeraDetectionVerbose = TRUE
			))

# Can extend and create with new("fastSingle", ...)
#' Inherits methods and attributes from fastq2otu class
#' Lists paramaters that are essential for analysing single-end data
#' @keyword internal
#' @export
setClass("fastSingle",
			slots = c(
				# === Required inputs for single-end data ===
				pathToData = "character"
			),
			prototype = list(
				pathToData = "NA_character_"
			),
			contains = "fastq2otu",
			validity = check_fastq2otu)

# Can extend and create with new("fastPaired", ...)
#' Inherits methods and attributes from fastq2otu class
#' Lists paramaters that are essential for analysing paired-end data
#' @keyword internal
#' @export
setClass("fastPaired",
			slots = c(
				# === Required inputs for paired-end data ===
				pathToData = "character",
				
				# === Merging Pairs ===
				doMergeSeqPairs = "logical",
				mergeSeqPairsTrimOverhang = "logical",
				mergeSeqPairsMinOverlap = "numeric",
				mergeSeqPairsMaxMismatch = "numeric", 
				mergeSeqPairsReturnRejects = "logical",
				mergeSeqPairsJustConcatenate = "logical", 
				mergeSeqPairsVerbose = "logical"
			),
			prototype = list(
				pathToData = "NA_character_",
											
				doMergeSeqPairs = FALSE,
				mergeSeqPairsTrimOverhang = FALSE,
				mergeSeqPairsMinOverlap = 12,
				mergeSeqPairsMaxMismatch = 0,
				mergeSeqPairsReturnRejects = FALSE,
				mergeSeqPairsJustConcatenate = FALSE,
				mergeSeqPairsVerbose = FALSE
			),
			contains = "fastq2otu",
			validity = check_fastq2otu)
			
			
# Can extend and create with new("fastFilter", ...)
#' FASTFILTER Custom Class
#'
#' Use to set parameters for DADA2 filterAndTrim function. Inputs for maxEE, truncQ, truncLen, and 
#' trimming functions (trimLeft, trimRight) are validated by check_filt_params function. 
#' 
#' Sub-class of fastSingle and fastPaired S4 datatypes.
#' @keyword internal
#' @export
setClass("fastFilter", 
			slots = c(
				# === Set parameters for filtering and trimming ===
				filtMaxEE = "numeric",
				filtTruncQ = "numeric",
				filtTruncLen = "numeric",
				filtTrimLeft = "numeric",
				filtTrimRight = "numeric",
				filtMultiThread = "logical",
				filtVerbose = "logical",
				filtMatchIDs = "logical",
				filtMinLen = "numeric"),
			prototype = list(
				filtMaxEE = 2.5,
				filtTruncQ = 0, 
				filtTruncLen = 0,
				filtTrimLeft = 0,
				filtTrimRight = 0,
				filtMultiThread = TRUE, 
				filtVerbose = TRUE,
				filtMatchIDs = FALSE, 
				filtMinLen = 50),	
			contains = c("fastSingle", "fastPaired"),
			validity = check_filt_params)


# Can extend and create with new("fastAssignTaxa", ...)
#' FASTASSIGNTAXA Custom Class
#'
#' Use to set parameters for DADA2 assignTaxonomy function. 
#' 
#' Sub-class of fastq2otu S4 datatype.
#' @keyword internal
#' @export
setClass("fastAssignTaxa", 
			slots = c(
				# === Set parameters for taxonomic assignement ==="
				taxDatabase = "character",
				assignTaxMinBootstrap = "numeric",
				assignTaxTryComplement = "logical",
				assignTaxOutputBootstraps = "logical",
				assignTaxLevels = "character",
				assignTaxMultiThread = "logical",
				assignTaxVerbose = "logical"),
			prototype = list(
				assignTaxMinBootstrap = 50,
				assignTaxTryComplement = FALSE,
				assignTaxOutputBootstraps = FALSE,
				assignTaxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
				assignTaxMultiThread = TRUE,
				assignTaxVerbose = TRUE),
			contains = c("fastSingle", "fastPaired"),
			validity = check_assign_tax
		)

# Can extend and create with new("fastPlotQuality", ...)
#' FASTPLOTQUALITY Custom Class
#'
#' Use to set parameters for DADA2's plotQualityProfile function.
#' 
#' Sub-class of fastq2otu S4 datatype.
#' @keyword internal
#' @export
setClass("fastPlotQuality", 
			slots = c(
				# === Generate a quality plot ===
				aggregateQual = "logical",
				qualN = "numeric"),
			prototype = list(
				aggregateQual = TRUE,
				qualN = 5e+05),
			contains = "fastq2otu"
		)

# Can extend and create with new("fastSeqDump", ...)
#' FASTSEQDUMP Custom Class
#'
#' Use to set parameters for retrieving FASTQ files from NCBI's SRA database
#' 
#' Sub-class of fastq2otu S4 datatype.
#' @keyword internal
#' @export
setClass("fastSeqDump", 
			slots = c(
				outDir = "character",
				pathToSampleURLs = "character",
				pathToSampleIDs = "character",
				pathToFastqDump = "character"),
			prototype = list(
				outDir = "NA_character_",
				pathToSampleURLs = "NA_character_",
				pathToSampleIDs = "NA_character_",
				pathToFastqDump = "NA_character_"),
			contains = "fastq2otu",
			validity = check_seq_dump
		)

# Can extend and create with new("fastPrimerTrim", ...)
#' FASTPRIMERTRIM Custom Class
#'
#' Use to set parameters for removing adapter sequences from amplicons using BBTools' bbduk.sh
#' script.
#' 
#' Sub-class of fastSingle and fastPaired S4 datatypes.
#' @keyword internal
#' @export
setClass("fastPrimerTrim",
			slots = c(
				# === Trim primers with bbduk.sh ===
				trimPrimers = "logical",
				listOfAdapters = "character",
				pathToRawFastq = "character",
				pathToNoPrimers = "character"),
			prototype = list(
				trimPrimers = FALSE,
				listOfAdapters = "NA_character_",
				pathToRawFastq = "NA_character_",
				pathToNoPrimers = "NA_character_"),
			contains = c("fastSingle", "fastPaired")
		)

# Can extend and create with new("fastReport", ...)
#' FASTREPORT Custom Class
#'
#' Use to set parameters for generating fastqc reports using FASTQCR.
#' Requires FASTQC to be installed onto system.
#' 
#' Sub-class of fastq2otu S4 datatype.
#' @keyword internal
#' @export
setClass("fastReport", 
			slots = c(
				# === Use Fastqc to generate data report ===
				runFastqc = "logical",
				installFastqc = "logical", 
				pathToFastqc = "character", 
				pathToFastqcResults = "character",
				fastqcThreads = "numeric",
				fastqcExperimentDescription = "character"),
			prototype = list(
				runFastqc = FALSE,
				installFastqc = FALSE, 
				pathToFastqc = "NA_character_", 
				pathToFastqcResults = "NA_character_",
				fastqcThreads = 4,
				fastqcExperimentDescription = "NA_character_"),
			contains = "fastq2otu"
		)

#===================
# Create class methods
#===================

#' Create a fastPaired object
#' 
#' @param inDir Required input. Path to FASTQ files. 
#' @param outDir Required input. Path to write output files.
#' @param mergeSeqs Required input. Default is TRUE to merge forward and reverse amplicons. FALSE prevents DADA2's mergePairs method from being executed.
#' @param trimOverhang Default is FALSE. Sets trimOverhang parameter for DADA2's mergePairs method. 
#' @param minOverlap Default is 12. Sets minOverlap parameter for DADA2's mergePairs function.
#' @param maxMismatch Default is 0. Sets maxMismatch parameter for DADA2's mergePairs method.  
#' @param returnRejects Default is FALSE. Sets returnRejects paramater for DADA2's mergePairs method.  
#' @param justConcatenate Default is FALSE. Sets returnRejects paramater for DADA2's mergePairs method. 
#' @param verbose Default is FALSE. Sets verbose parameter for DADA2 methods. 
#' @param multithread Default is FALSE. Sets multithread parameter for DADA2 methods. 
#' @param prefix Default is "myproject". Project-specific label that will be included in the names of all files and directories created.
#' @param isPaired Default is TRUE for fastPaired object. If FALSE, a fastSingle object will be created instead. 
#' @param derepN Default is 1e+06. Sets the sampling parameter for DADA2's derepFastq method.
#' @param getErrPDF Default is FALSE. If TRUE a learn errors plot from DADA's plotErrors method is produced. 
#' @param errN Default is 1e+08. Sets the nbases parameter for DADA2's learnErrors method.  
#' @param dadaBandSize Default is 16. Sets the BAND_SIZE parameter for DADA2's dada method. 
#' @param dadaOmegaA Default is 1e-40. Sets the OMEGA_A paramter for DADA2's dada method.
#' @param getChimeraTable Default is FALSE. When TRUE, creates table showing all sequences labled as chimeric by DADA2 (uses DADA2's isBimeraDevovoTable function as opposed to removeBimeraDenovo). 
#' @param minSampleFraction Default is 0.9. Sets minDetection parameter for DADA2's removeBimeraDenovo and isBimera functions.
#' @param ignoreNegatives Default is 1. Sets ignoreNegatives parameter for DADA2's removeBimeraDenovo and isBimera functions. 
#' @param minFoldParentOverAbundance Default is 1.5. Sets minfoldParentOverabundance parameter for DADA2's removeBimeraDenovo and isBimera functions. 
#' @param detectionAbundance Default is 2. Sets detectionParentAbundace parameter for DADA2's removeBimeraDenovo and isBimera functions. 
#' @param allowOneOff Default is FALSE. Sets allowOneOff parameter for DADA2's removeBimeraDenovo and isBimera functions. 
#' @param maxShift Default is 16. Sets maxShift parameter for DADA2's removeBimeraDenovo and isBimera functions. 
#' @param downloadSeqs. Default is FALSE. If TRUE, a fastSeqDump object will be created when runPipeline method is executed.
#' @param trimAdapters. Default is FALSE. If TRUE, a fastPrimerTrim object will be created when runPipeline method is executed. 
#' @param generateReport Default is FALSE. If TRUE, a fastReport object will be created when runPipeline is executed. 
#' 
#  @return S4 object of type fastPaired
#' @export
setFastPaired <- function(inDir, outDir, mergeSeqs = TRUE, trimOverhang = FALSE, minOverlap = 12,
					maxMismatch = 0, returnRejects = FALSE, justConcatenate = FALSE, verbose = FALSE, 
					prefix = "myproject", isPaired = TRUE, derepN = 1e+06, getErrPDF = FALSE, 
					errN = 1e+08, multithread = FALSE, dadaBandSize = 16, dadaOmegaA = 1e-40, 
					getChimeraTable = FALSE, minSampleFraction = 0.9, ignoreNegatives = 1, 
					minFoldParentOverAbundance = 1.5, detectionAbundance = 2, allowOneOff = FALSE, maxShift = 16,
					downloadSeqs = FALSE, trimAdapters = FALSE, generateReport = FALSE) {
			
			temp <- methods::new("fastPaired", 
						pathToData = inDir, 
						outDir = outDir,
						projectPrefix = prefix, 
						isPaired = isPaired, 
												
						downloadSeqs = downloadSeqs, 
						trimAdapters = trimAdapters,
						generateReport = generateReport,
						
						doMergeSeqPairs = mergeSeqs,
						mergeSeqPairsTrimOverhang = trimOverhang,
						mergeSeqPairsMinOverlap = minOverlap,
						mergeSeqPairsMaxMismatch = maxMismatch, 
						mergeSeqPairsReturnRejects = returnRejects,
						mergeSeqPairsJustConcatenate = justConcatenate, 
						mergeSeqPairsVerbose = verbose,
												
						derepVerbose = verbose,
						derepN = derepN,

						errN = errN,
						saveErrorsPlotPDF = getErrPDF, 
						errMultithread = multithread,

						dadaBandSize = dadaBandSize,
						dadaOmegaA = dadaOmegaA,

						createChimeraDetectionTable = getChimeraTable,
						chimeraDetectionMinSampleFraction = minSampleFraction,
						chimeraDetectionIgnoreNegatives = ignoreNegatives,
						chimeraDetectionMinFoldParentOverabundance = minFoldParentOverAbundance,
						chimeraDetectionParentAbundance = detectionAbundance,
						chimeraDetectionAllowOneOff = allowOneOff,
						chimeraDetectionMaxShift = maxShift)
			
			return(temp)
						
	
}

#' Create a fastSingle object
#' 
#' @param inDir Required input. Path to FASTQ files. 
#' @param outDir Required input. Path to write output files.
#' @param verbose Default is FALSE. Sets verbose parameter for DADA2 methods. 
#' @param multithread Default is FALSE. Sets multithread parameter for DADA2 methods. 
#' @param prefix Default is "myproject". Project-specific label that will be included in the names of all files and directories created.
#' @param isPaired Default is TRUE for fastPaired object. If FALSE, a fastSingle object will be created instead. 
#' @param derepN Default is 1e+06. Sets the sampling parameter for DADA2's derepFastq method.
#' @param errN Default is 1e+08. Sets the nbases parameter for DADA2's learnErrors method. 
#' @param getErrPDF Default is FALSE. If TRUE a learn errors plot from DADA's plotErrors method is produced.  
#' @param dadaBandSize Default is 16. Sets the BAND_SIZE parameter for DADA2's dada method. 
#' @param dadaOmegaA Default is 1e-40. Sets the OMEGA_A paramter for DADA2's dada method.
#' @param getChimeraTable Default is FALSE. When TRUE, creates table showing all sequences labled as chimeric by DADA2 (uses DADA2's isBimeraDevovoTable function as opposed to removeBimeraDenovo). 
#' @param minSampleFraction Default is 0.9. Sets minDetection parameter for DADA2's removeBimeraDenovo and isBimera functions.
#' @param ignoreNegatives Default is 1. Sets ignoreNegatives parameter for DADA2's removeBimeraDenovo and isBimera functions. 
#' @param minFoldParentOverAbundance Default is 1.5. Sets minfoldParentOverabundance parameter for DADA2's removeBimeraDenovo and isBimera functions. 
#' @param detectionAbundance Default is 2. Sets detectionParentAbundace parameter for DADA2's removeBimeraDenovo and isBimera functions. 
#' @param allowOneOff Default is FALSE. Sets allowOneOff parameter for DADA2's removeBimeraDenovo and isBimera functions. 
#' @param maxShift Default is 16. Sets maxShift parameter for DADA2's removeBimeraDenovo and isBimera functions. 
#' @param downloadSeqs. Default is FALSE. If TRUE, a fastSeqDump object will be created when runPipeline method is executed.
#' @param trimAdapters. Default is FALSE. If TRUE, a fastPrimerTrim object will be created when runPipeline method is executed. 
#' @param generateReport Default is FALSE. If TRUE, a fastReport object will be created when runPipeline is executed. 
#' 
#' @return S4 object of type fastSingle
#' @export
setFastSingle <- function(inDir, outDir, verbose = FALSE, 
					prefix = "myproject", isPaired = FALSE, derepN = 1e+06, getErrPDF = FALSE, 
					errN = 1e+08, multithread = FALSE, dadaBandSize = 16, dadaOmegaA = 1e-40, 
					getChimeraTable = FALSE, minSampleFraction = 0.9, ignoreNegatives = 1, 
					minFoldParentOverAbundance = 1.5, detectionAbundance = 2, allowOneOff = FALSE, maxShift = 16,
					downloadSeqs = FALSE, trimAdapters = FALSE, generateReport = FALSE) {
			
				temp <- methods::new("fastSingle", 
						pathToData = inDir, 
						outDir = outDir,
						projectPrefix = prefix, 
						isPaired = isPaired, 
												
						downloadSeqs = downloadSeqs, 
						trimAdapters = trimAdapters,
						generateReport = generateReport,
																		
						derepVerbose = verbose,
						derepN = derepN,

						errN = errN,
						saveErrorsPlotPDF = getErrPDF, 
						errMultithread = multithread,

						dadaBandSize = dadaBandSize,
						dadaOmegaA = dadaOmegaA,

						createChimeraDetectionTable = getChimeraTable,
						chimeraDetectionMinSampleFraction = minSampleFraction,
						chimeraDetectionIgnoreNegatives = ignoreNegatives,
						chimeraDetectionMinFoldParentOverabundance = minFoldParentOverAbundance,
						chimeraDetectionParentAbundance = detectionAbundance,
						chimeraDetectionAllowOneOff = allowOneOff,
						chimeraDetectionMaxShift = maxShift)
			
			return(temp)		
}


#' Create a fastReport object
#' Users have the option of creating a detailed report of their data using FASTQC
#' 
#' @param inDir Required input. Path to directory containing FASTQ files. 
#' @param outDir Required input. Path to write output files.
#' @param adapterList Required input. List of adapter sequences to remove from files. 
#' 
#' TODO: Add Example 
#' 
#  @return S4 object of type fastReport
#' @export
setFastReport <- function(inDir, outDir, fastqcPath = NA_character_, installFastqc = FALSE,
					numThreads = NA, description = NA_character_) {
	temp <- new("fastReport",
				installFastqc = installFastqc, 
				pathToFastqc = fastqcPath, 
				pathToFastqcResults = outDir,
				fastqcThreads = numThread,
				fastqcExperimentDescription = description)
	return(temp)			
}

#' Create a fastPrimerTrim object
#' Users have the option of trimming adapters from their data using bbduk.sh 
#' 
#' @param inDir Required input. Path to directory containing FASTQ files. 
#' @param outDir Required input. Path to write output files.
#' @param adapterList Required input. List of adapter sequences to remove from files. 
#' 
#' TODO: Add Example 
#' 
#  @return S4 object of type fastPrimerTrim
#' @export
setFastPrimerTrim <- function(inDir, outDir, adapterList) {
	temp <- new("fastPrimerTrim", 
				listOfAdapters = adapterList,
				pathToRawFastq = inDir,
				pathToNoPrimers = outDir)
	
	return(temp)
}

#' Create a fastSeqDump object
#' Users have the option of downloading SRA data using fastq-dump or FTP. 
#' 
#' @param sampleURLs OPTIONAL input. Path to text file containing FTP URLs generated from SRA Explorer. If this parameter is provided, then fastqDumpPath is ignored.
#' @param outDir Required input. Path to write output files.
#' @param sampleList Required input if fastqDumpPath is specified. New-line delimited text file containin SRA accession IDs. 
#' @param fastqDumpPath Set path to fastq-dump script. If provided, fastq-dump is used in place of wget/curl.
#' 
#' TODO: Add Example 
#' 
#  @return S4 object of type fastSeqDump
#' @export
setFastSeqDump <- function(sampleURLs, outDir, sampleList, fastqDumpPath = NA_character_) {
	temp <- new("fastSeqDump", 
				outDir = outDir, 
				pathToSampleURLs = sampleURLs,
				pathToSampleIDs = sampleList,
				pathToFastqDump = fastqDumpPath)
				
	return(temp)
}

#' Create a fastFilter object
#' 
#' @param inDir Required input. Path to FASTQ files. 
#' @param outDir Required input. Path to write output files.
#' @param prefix Default is "myproject". Project-specific label that will be included in the names of all files and directories created.
#' 
#' @param verbose Default is FALSE. Sets verbose parameter for DADA2 methods. 
#' @param maxEE Default is 2.5. Sets maxEE parameter for DADA2's filterAndTrim method. Two numeric inputs are expected for paired-end data.
#' @param truncQ Default is 0. Sets truncQ parameter for DADA2's filterAndTrim method. Two numeric inputs are expected for paired-end data.
#' @param truncLen Default is 0. Sets truncLen parameter for DADA2 filterAndTrim function. Two numeric inputs are expected for paired-end data.
#' @param trimLeft Default is 0. Sets trimLeft parameter for DADA2 filterAndTrim function. Two numeric inputs are expected for paired-end data.
#' @param trimRight Default is 0. Sets trimRight parameter for DADA2 filterAndTrim function. Two numeric inputs are expected for paired-end data.
#' @param matchIDs Default is FALSE. Sets matchIDs parameter for DADA2 filterAndTrim function. Set TRUE when handling paired-end data with complementary file names for forward and reverse reads.
#' @param minLen Default is 50. Sets minLen parameter for DADA2 filterAndTrim function. Two numeric inputs are expected for paired-end data.
#' @param multithread Default is FALSE. Sets multithread parameter for DADA2 methods. 
#' @param isPaired Default is FALSE. If TRUE, two values must be provided for truncQ, truncLen, trimLeft and trimRight, and minLen (see example below).
#' 
#' TODO: Add Example 
#' 
#  @return S4 object of type fastFilter
#' @export
setFastFilter <- function(inDir, outDir, prefix = "myproject", 
					maxEE = 2.5, truncQ = 0, truncLen = 0, trimLeft = 0,
					trimRight = 0, matchIDs = FALSE, minLen = 50,
					isPaired = FALSE, verbose = FALSE, multithread = FALSE) {
				
				
			temp <- methods::new("fastFilter", 
					pathToData = inDir, 
					outDir = outDir,
					isPaired = isPaired, 
					projectPrefix = prefix, 
									
					filtMaxEE = maxEE,
					filtTruncQ = truncQ,
					filtTruncLen = truncLen,
					filtTrimLeft = trimLeft,
					filtTrimRight = trimRight,
					filtMatchIDs = matchIDs,
					filtVerbose = verbose,
					filtMultiThread = multithread,
					filtMinLen = minLen)
			
			return(temp)
}

#' Create a fastAssignTaxa object
#' 
#' @param refDatabase Required input. Path to fasta-formatted database containing reference sequences used during taxonomic assignment. 
#' @param prefix Default is "myproject". Project-specific label that will be included in the names of all files and directories created.
#' 
#' @param minBootstrap Default is 50. Sets min-bootstrap value for DADA2 assignTaxonomy method.
#' @param tryComplement Default is FALSE. Sets value for tryComplement parameter for DADA2's assignTaxonomy method. 
#' @param showBootstraps Default is FALSE. Sets value for showBootstraps parameter for DADA2's assignTaxonomy method.
#' @param taxLevels Provide character list of taxonomic levels to include in taxonomy assignment. The default OTU table generated from this function includes assignments from kingdom to species.
#'
#' @param verbose Default is FALSE. Sets verbose parameter for DADA2 methods. 
#' @param multithread Default is FALSE. Sets multithread parameter for DADA2 methods. 
#' 
#' TODO: Add Example 
#' 
#  @return S4 object of type fastAssignTaxa
#' @export
setFastAssignTaxa <- function(refDatabase, prefix = "myproject",
					minBootstrap = 50, tryComplement = FALSE, showBootstraps = FALSE,
					taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
					verbose = FALSE, multithread = FALSE
					) {
					
			temp <- methods::new("fastAssignTaxa", 
				projectPrefix = prefix, 
								
				taxDatabase = refDatabase,
				assignTaxMinBootstrap = minBootstrap,
				assignTaxTryComplement = tryComplement,
				assignTaxOutputBootstraps = showBootstraps,
				assignTaxLevels = taxLevels,
				assignTaxMultiThread = multithread,
				assignTaxVerbose = verbose)
			
			return(temp)

}

#' Create a fastPlotQuality object
#' 
#' @param aggregate Default if FALSE. Sets aggregate parameter for DADA's plotQualityProfile method. 
#' @param N Default is 5e+05. Sets N parameter for DADA's plotQualityProfile method.
#' 
#' TODO: Add Example 
#' 
#  @return S4 object of type fastPlotQuality
#' @export
setfastPlotQuality <- function(aggregate = FALSE, N = 5e+05) {
			temp <- methods::new("fastPlotQuality", 
				aggregateQual = aggregate,
				qualN = N)
}





