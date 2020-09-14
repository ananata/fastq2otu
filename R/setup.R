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
  if (object@checkValidity) {
  	errors <- character()
    length_output <- length(object@outDir)
  	if (length_output == 0) {
  		msg <- paste("No input provided for outDir. A value is required")
  		errors <- c(errors, msg)
  	}
    else if (!is.na(object@outDir) & !dir.exists(object@outDir)) {
      msg <- paste("Provided path for outDir: ", object@outDir, " is invalid and does not exist", sep = "")
      errors <- c(errors, msg)
    }
  	if (length(errors) == 0) TRUE else errors
  }
}

#' Check for valid inputs for fastPaired
#'
#' This function checks to see if the isPaired parameter reflects the object the user is attempting to create.
#' Generates an additional level of validity that insures uses is using package correctly.
#'
#' @param object fastPaired-type S4 object
#' @return TRUE (boolean) or error message detailing source(s) of error.
#' @keyword internal
#' @export
check_fastSingle <- function(object) {
  if (object@checkValidity) {
    errors <- character()
    if (object@isPaired) {
      msg <- "fastSingle object cannot be paired. Please create a fastPaired object instead"
      errors <- c(errors, msg)
    }
    if (object@pathToData == "" | is.na(object@pathToData)) { # checks for empty strings
      msg <- paste("No path provided for pathToData.")
      errors <- c(errors, msg)
    }
    else if (!is.na(object@pathToData) & !dir.exists(object@pathToData)) {
      msg <- paste("Provided path for pathToData: ", object@pathToData, " is invalid and does not exist", sep = "")
      errors <- c(errors, msg)
    }
    
    if (length(errors) == 0) TRUE else errors
  }
}

#' Check for valid inputs for fastPaired
#'
#' This function checks to see if the isPaired parameter reflects the object the user is attempting to create.
#' Generates an additional level of validity that insures uses is using package correctly.
#'
#' @param object fastPaired-type S4 object
#' @return TRUE (boolean) or error message detailing source(s) of error.
#' @keyword internal
#' @export
check_fastPaired <- function(object) {
  if (object@checkValidity) {
    errors <- character()
    if (!object@isPaired) {
      msg <- "fastPaired object must use paired-end data. Please create a fastSingle object instead"
      errors <- c(errors, msg)
    }
    if (object@pathToData == "" | is.na(object@pathToData)) { # checks for empty strings
      msg <- paste("No path provided for pathToData.")
      errors <- c(errors, msg)
    }
    else if (!is.na(object@pathToData) & !dir.exists(object@pathToData)) {
      msg <- paste("Provided path for pathToData: ", object@pathToData, " is invalid and does not exist", sep = "")
      errors <- c(errors, msg)
    }
    if (length(errors) == 0) TRUE else errors
  }
}


#' Check for valid inputs for fastPrimerTrim
#'
#' This function checks to see if "pathToBBduk", "listOfAdapters", "pathToRawFastq" and "pathToNoPrimers" parameters are present and valid. Missing or invalid inputs will cause an
#' error message to be displayed.
#'
#' @param object fastPrimerTrim-type S4 object
#' @return TRUE (boolean) or error message detailing source(s) of error.
#' @keyword internal
#' @export
check_primer_trim <- function(object) {
	errors <- character()
	if (is.na(object@listOfAdapters) | object@listOfAdapters == "") {
		msg <- "No valid input provided for listOfAdapters. Please update your config file and try again"
		errors <- c(errors, msg)
	} else if (!file.exists(object@listOfAdapters)) {
		msg <- paste0("The path provided for listOfAdapters: ", object@listOfAdapters, " is invalid and does not exist")
		errors <- c(errors, msg)
	}
	if (is.na(object@pathToRawFastq) | object@pathToRawFastq == "") {
		msg <- "No valid input provided for pathToRawFastq. Please update your config file and try again. "
		errors <- c(errors, msg)
	} else if (!dir.exists(object@pathToRawFastq)) {
		msg <- paste0("The path provided for pathToRawFastq: ", object@pathToRawFastq, " is invalid and does not exist")
		errors <- c(errors, msg)
	}
	if (is.na(object@pathToNoPrimers) | object@pathToNoPrimers == "") {
		msg <- "No valid input provided for pathToNoPrimers. Please update your config file and try again. "
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
		if (length(object@filtTruncQ) != 2) { 
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
		if (length(object@filtTruncQ) != 1) { 
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
	if (is.null(object@pathToSampleURLs) & !is.null(object@pathToSampleIDs)) {
	      msg1 <- "Unable to find pathToSampleURLs and pathToSampleIDs parameters"
	}
	
	if (!file.exists(object@pathToSampleURLs) & !file.exists(object@pathToSampleIDs)) {
				msg1 <- paste0("Value for pathToSampleURLs: ", object@pathToSampleURLs, " could not be found. Please enter valid file path.")
                        errors <- c(errors, msg1)
				msg2 <- paste0("Value for pathToSampleIDs: ", object@pathToSampleIDs, " could not be found. Please enter valid file path.")
                        errors <- c(errors, msg2)
	} else {
	  if (file.exists(object@pathToSampleIDs) & object@installFastqDump) {
      if (is.null(object@pathToFastqDump)) {
        msg <- "pathToFastqDump parameter must be present and valid when installFastqDump parameter is TRUE"
        errors <- c(errors, msg)
      }
      if (!is.null(object@pathToFastqDump) & !file.exists(object@pathToFastqDump)) {
        msg <- "Invalid input for pathToFastqDump. Please provide a valid path to fastq-dump"
        errors <- c(errors, msg)
      }
	  }
	} 
	if (length(errors) == 0) TRUE else errors
}

#' Check for valid inputs for fastAssignTaxa
#'
#' This function checks to see if fastAssignTaxa parameters are present and valid. Missing or invalid inputs will cause an
#' error message to be displayed.
#'
#' @param object fastAssignTaxa-type S4 object
#' @return TRUE (boolean) or error message detailing source(s) of error.
#' @keyword internal
#' @export
check_assign_tax <- function(object) {
	errors <- c()
	if (!file.exists(object@taxDatabase)) {
		if (is.na(object@taxDatabase) | length(object@taxDatabase) == 0 ) {
			msg <- ("No path to reference database provided.")
			errors <- c(errors, msg)
		} else {
			msg <- paste("Path to database: ", object@taxDatabase, " is invalid and does not exist", sep = "")
			errors <- c(errors, msg)
		}
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
					checkValidity = "logical",
					
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
				outDir = NA_character_,
				isPaired = FALSE, 
				checkValidity = FALSE, 

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
			),
		validity = check_fastq2otu)

# Can extend and create with new("fastSingle", ...)
#' Inherits methods and attributes from fastq2otu class
#' Lists paramaters that are essential for analysing single-end data
#' @keyword internal
#' @export
setClass("fastSingle",
			slots = c(
				# === Required inputs for single-end data ===
				pathToData = "character",
				filePattern = "character"
			),
			prototype = list(
				pathToData = NA_character_,
				filePattern = "*.fastq(.gz)?$"
			),
			contains = "fastq2otu",
			validity = check_fastSingle)

# Can extend and create with new("fastPaired", ...)
#' Inherits methods and attributes from fastq2otu class
#' Lists paramaters that are essential for analysing paired-end data
#' @keyword internal
#' @export
setClass("fastPaired",
			slots = c(
				# === Required inputs for paired-end data ===
				pathToData = "character",
				filePattern = "character",
				
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
				pathToData = NA_character_,
				filePattern = ".*_[1,2].fastq(.gz)?",											
				doMergeSeqPairs = FALSE,
				mergeSeqPairsTrimOverhang = FALSE,
				mergeSeqPairsMinOverlap = 12,
				mergeSeqPairsMaxMismatch = 0,
				mergeSeqPairsReturnRejects = FALSE,
				mergeSeqPairsJustConcatenate = FALSE,
				mergeSeqPairsVerbose = FALSE
			),
			contains = "fastq2otu",
			validity = check_fastPaired)
			
			
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
				# === Set parameters for taxonomic assignment ==="
				taxDatabase = "character",
				assignTaxMinBootstrap = "numeric",
				assignTaxTryComplement = "logical",
				assignTaxOutputBootstraps = "logical",
				assignTaxLevels = "character",
				assignTaxMultiThread = "logical",
				assignTaxVerbose = "logical"),
			prototype = list(
				taxDatabase = NA_character_,
				assignTaxMinBootstrap = 50,
				assignTaxTryComplement = FALSE,
				assignTaxOutputBootstraps = FALSE,
				assignTaxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
				assignTaxMultiThread = TRUE,
				assignTaxVerbose = TRUE),
			contains = "fastq2otu",
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
				pathToFastqDump = "character",
				installFastqDump = "logical", 
				useDump = "logical"),
			prototype = list(
				outDir = NA_character_,
				pathToSampleURLs = NA_character_,
				pathToSampleIDs = NA_character_,
				pathToFastqDump = NA_character_,
				installFastqDump = FALSE,
				useDump = FALSE),
			contains = c("fastSingle", "fastPaired"),
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
				pathToNoPrimers = "character",
				maxMismatch = "numeric", 
				allowIndels = "logical",
				trimForwardReads = "logical",
				trimReverseReads = "logical",
				orientReads = "logical",
				compressFiles = "logical"),
			prototype = list(
				trimPrimers = FALSE,
				listOfAdapters = NA_character_,
				pathToRawFastq = NA_character_,
				pathToNoPrimers = NA_character_,
				maxMismatch = 2,
				allowIndels = FALSE,
				trimForwardReads = FALSE,
				trimReverseReads = FALSE,
				orientReads = FALSE, 
				compressFiles = FALSE),
			contains = "fastq2otu",
			validity = check_primer_trim
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
				pathToFastqc = NA_character_, 
				pathToFastqcResults = NA_character_,
				fastqcThreads = 4,
				fastqcExperimentDescription = NA_character_),
			contains = "fastq2otu"
		)

#===================
# Create class methods
#===================

#' Create a fastPaired object
#' 
#' @param inDir Required input. Path to FASTQ files. 
#' @param outDir Required input. Path to write output files.
#' @param validate Default is TRUE. Execute all validation methods.
#' @param pattern Regex pattern used to select fastq files from directory.
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
setFastPaired <- function(inDir, outDir, validate = TRUE, pattern = ".*_[1,2].fastq(.gz)?", mergeSeqs = TRUE, trimOverhang = FALSE, minOverlap = 12,
					maxMismatch = 0, returnRejects = FALSE, justConcatenate = FALSE, verbose = FALSE, 
					prefix = "myproject", isPaired = TRUE, derepN = 1e+06, getErrPDF = FALSE, 
					errN = 1e+08, multithread = FALSE, dadaBandSize = 16, dadaOmegaA = 1e-40, 
					getChimeraTable = FALSE, minSampleFraction = 0.9, ignoreNegatives = 1, 
					minFoldParentOverAbundance = 1.5, detectionAbundance = 2, allowOneOff = FALSE, maxShift = 16,
					downloadSeqs = FALSE, trimAdapters = FALSE, generateReport = FALSE) {
			
			temp <- methods::new("fastPaired", 
						pathToData = inDir, 
						outDir = outDir,
						filePattern = pattern, 
						projectPrefix = prefix, 
						isPaired = isPaired, 
						checkValidity = validate,
												
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
#' @param validate Default is TRUE. Executes normal validation steps
#' @param pattern Regex pattern used to select fastq files from directory.
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
setFastSingle <- function(inDir, outDir, validate = TRUE, verbose = FALSE, pattern = ".*.fastq(.gz)?",
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
						checkValidity = validate, 
						filePattern = pattern, 
												
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
#' @param validate Default is FALSE. When TRUE, all validation functions are executed.
#' @param inDir Required input. Path to directory containing FASTQ files. 
#' @param outDir Required input. Path to write output files.
#' @param adapterList Required input. List of adapter sequences to remove from files. 
#' 
#' TODO: Add Example 
#' 
#  @return S4 object of type fastReport
#' @export
setFastReport <- function(inDir, outDir, validate = FALSE, fastqcPath = NA_character_, installFastqc = FALSE,
					numThreads = 4, description = NA_character_) {
	temp <- new("fastReport",
	      checkValidity = validate,
				installFastqc = installFastqc, 
				pathToFastqc = fastqcPath, 
				pathToFastqcResults = outDir,
				fastqcThreads = numThreads,
				fastqcExperimentDescription = description)
	return(temp)			
}

#' Create a fastPrimerTrim object
#' Users have the option of trimming adapters from their data using bbduk.sh 
#' 
#' @param inDir Required input. Path to directory containing FASTQ files. 
#' @param outDir Required input. Path to write output files.
#' @param adapterList Required input. List of adapter sequences to remove from files. 
#' @param validate Default is FALSE. When TRUE, all validation method are executed.
#' @param paired Default is FALSE. When TRUE, two adapter sequences will be expected in file listed in "listOfAdapters"
#' @param maxMismatch Default is 2. Sets the number of bases to mismatch when aligning forward / reverse adapters to sequences.
#' @param allowIndels Default is false. Allows gaps / insertions when aligning forward / reverse adapters to sequences. 
#' @param trimF Default is FALSE. Trim forward reads.
#' @param trimR Default is FALSE. Trim reverse reads.
#' @param orient Default is FALSE. Orient forward and reverse reads into the same direction
#' @param compress Default is FALSE. When TRUE, all output files are gzipped.
#' TODO: Add Example 
#' 
#  @return S4 object of type fastPrimerTrim
#' @export
setFastPrimerTrim <- function(inDir, outDir, adapterList, paired = FALSE, validate = FALSE,
			maxMismatch = 2, allowIndels = FALSE, trimF = FALSE, trimR = FALSE, 
			orient = FALSE, compress = FALSE) {
	temp <- new("fastPrimerTrim", 
	      checkValidity = validate,
				listOfAdapters = adapterList,
				pathToRawFastq = inDir,
				pathToNoPrimers = outDir,
				isPaired = paired, 
				maxMismatch = maxMismatch,
                                allowIndels = allowIndels,
                                trimForwardReads = trimF,
                                trimReverseReads = trimR,
                                orientReads = orient,
                                compressFiles = compress)
	
	return(temp)
}

 #' Create a fastSeqDump object
#' Users have the option of downloading SRA data using fastq-dump or FTP. 
#' 
#' @param sampleURLs OPTIONAL input. Path to text file containing FTP URLs generated from SRA Explorer. If this parameter is provided, then fastqDumpPath is ignored.
#' @param outDir Required input. Path to write output files.
#' @param validate Default is FALSE. When TRUE, all validation methods are executed
#' @param sampleList Required input if fastqDumpPath is specified. New-line delimited text file containin SRA accession IDs. 
#' @param useFastqDump If TRUE, Fastq-dump will be used in place of wget. Requires fastqDumpPath to be supplied. If not present, fastq-dump present in system.files will be used.
#' @param fastqDumpPath Set path to fastq-dump script. If provided, fastq-dump is used in place of wget/curl.
#' 
#' TODO: Add Example 
#' 
#  @return S4 object of type fastSeqDump
#' @export
setFastSeqDump <- function(sampleURLs, outDir, validate = FALSE, sampleList, useFastqDump = FALSE, fastqDumpPath = NA_character_, install_dump = FALSE) {
	temp <- new("fastSeqDump", 
				outDir = outDir, 
				checkValidity = validate,
				pathToSampleURLs = sampleURLs,
				pathToSampleIDs = sampleList,
				useDump = useFastqDump, 
				installFastqDump = install_dump,
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
#' @param validate Default is FALSE. When TRUE all validation methods are executed. 
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
setFastFilter <- function(inDir, outDir, prefix = "myproject", validate = FALSE,  
					maxEE = 2.5, truncQ = 0, truncLen = 0, trimLeft = 0,
					trimRight = 0, matchIDs = FALSE, minLen = 50,
					isPaired = FALSE, verbose = FALSE, multithread = FALSE) {
				
				
			temp <- methods::new("fastFilter", 
					pathToData = inDir, 
					outDir = outDir,
					isPaired = isPaired, 
					projectPrefix = prefix, 
					checkValidity = validate,
									
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
#' @param validate Default is FALSE. When TRUE all validation methods are executed. 
#' @param verbose Default is FALSE. Sets verbose parameter for DADA2 methods. 
#' @param multithread Default is FALSE. Sets multithread parameter for DADA2 methods. 
#' 
#' TODO: Add Example 
#' 
#  @return S4 object of type fastAssignTaxa
#' @export
setFastAssignTaxa <- function(refDatabase, prefix = "myproject", validate = FALSE,
					minBootstrap = 50, tryComplement = FALSE, showBootstraps = FALSE,
					taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
					verbose = FALSE, multithread = FALSE
					) {
					
			temp <- methods::new("fastAssignTaxa", 
				taxDatabase = refDatabase,
				projectPrefix = prefix, 
				checkValidity = validate, 
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
#' @param validate Default is FALSE, when TRUE all validation methods are executed. 
#' TODO: Add Example 
#' 
#  @return S4 object of type fastPlotQuality
#' @export
setfastPlotQuality <- function(aggregate = FALSE, N = 5e+05, validate = FALSE) {
			temp <- methods::new("fastPlotQuality", 
				aggregateQual = aggregate,
				checkValidity = validate, 
				qualN = N)
}





