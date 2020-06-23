#' Check for valid inputs for fastq2otu
#'
#' This function checks to see if user inputs are present and valid. Missing or invalid inputs will cause an
#' error message to be displayed.
#'
#' @param object fastq2otu-type S4 object
#' @return TRUE (boolean) or error message detailing source(s) of error.
#' @export
check_fastq2otu <- function(object) {
 	errors <- character()
  	if (!file.exists(object@taxDatabase)) {
    		msg <- paste("Provided input: ", object@taxDatabase, " is invalid and does not exist", sep = "")
    		errors <- c(errors, msg)
  	}
  	if (!dir.exists(object@pathToData)) {
    		msg <- paste("Provided input: ", object@pathToData, " is invalid and does not exist", sep = "")
   		errors <- c(errors, msg)
  	}
	length_prefix <- length(object@projectPrefix)
	if (length_prefix == 0) {
		msg <- paste("No input provided for projectPrefix. A value is required")
		errors <- c(errors, msg)
	}
        length_output <- length(object@outDir)
	if (length_output == 0) {
		msg <- paste("No input provided for outDir. A value is required")
		errors <- c(errors, msg)
	}
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
	if (length(errors) == 0) TRUE else errors
}


# Can extend and create with new("fastq2otu", ...)
#' Create a custom class for FASTQ2OTU 
#'
#' Used to convert parameters stored in YML config file to accessible attributes associated with
#' an S4-type data class. 
#' The object consists of elements or parameters that are useful for both single-end and paired-end datasets. 
#' It also list elements that are essential for executing Fastq-dump (NCBI SRA-toolkit), bbduk.sh (BBTools), and FASTQC
#' should a user desire. 
#' @keyword internal
#' @export
setClass("fastq2otu",
                slots = c(
			# === Required Inputs ===
			taxDatabase = "character",
			projectPrefix = "character",
			outDir = "character",
			isPaired = "logical", 

			# === Use Fastq-dump to download SRA data ===
			pathToSampleIDs = "character",
			runFastqDump = "logical",
			
			# === Trim primers with bbduk.sh ===
			trimPrimers = "logical"
			listOfAdapters = "character",
			pathToRawFastq = "character",
			pathToNoPrimers = "character",
			
			# === Use Fastqc to generate data report ===
			runFastqc = "logical",
			installFastqc = "logical", 
			pathToFastqc = "character", 
			pathToFastqcResults = "character",
			fastqcThreads = "numeric",
			fastqcExperimentDescription = "character",

			# === Create summary table that displays changes to total read counts ===
			finalSummaryTable = "character",
                      
			# === Generate a quality plot ===
			plotQuality = "logical",
			qualityPlotPDF = "character",

			# === Dereplicate reads to keep only unique sequences ===
			derepVerbose = "logical",
			derepN = "numeric",

			# === Detect and learn error patterns in sequences ===
			errN = "numeric",
			errMultithread = "numeric",
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
			chimeraDetectionVerbose = "logical",

			# === Set parameters for taxonomic assignement ==="
			assignTaxMinBootstrap = "numeric",
			assignTaxTryComplement = "logical",
			assignTaxOutputBootstraps = "logical",
			assignTaxLevels = "character",
			assignTaxMultiThread = "logical",
			assignTaxVerbose = "logical",
	
			# === Merge sample tables across to allow for cross sample comparisons ===
			mergeSamples = "logical",
			finalMergedTable = "character"
			),
		prototype = list(
			taxDatabase = "NA_character_",
			projectPrefix = "myproject",
			outDir = "NA_character_",
			isPaired = FALSE, 

			pathToSampleIDs = "NA_character_",
			runFastqDump = FALSE,

			trimPrimers = FALSE
			listOfAdapters = "NA_character_",
			pathToRawFastq = "NA_character_",
			pathToNoPrimers = "NA_character_",

			runFastqc = FALSE,
			installFastqc = FALSE, 
			pathToFastqc = "NA_character_", 
			pathToFastqcResults = "NA_character_",
			fastqcThreads = 4,
			fastqcExperimentDescription = "NA_character_",

			finalSummaryTable = "final_summary_table.txt",

			plotQuality = TRUE,
			qualityPlotPDF = "quality_plots.pdf",

			derepVerbose = TRUE,
			derepN = 1e+06,

			errN = 1e+08,
			errMultithread = 4,
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
			chimeraDetectionVerbose = TRUE,

			assignTaxMinBootstrap = 50,
			assignTaxTryComplement = FALSE,
			assignTaxOutputBootstraps = FALSE,
			assignTaxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
			assignTaxMultiThread = TRUE,
			assignTaxVerbose = TRUE,

			mergeSamples = TRUE,
			finalMergedTable = "final_merged_table.csv"
			),
		validity = check_fastq2otu)

# Can extend and create with new("fastq2otu_single", ...)
#' Inherits methods and attributes from fastq2otu class
#' Lists paramaters that are essential for analysing single-end data
#' @export
setClass("fastq2otu_single",
			slots = c(
				# === Required inputs for single-end data ===
				pathToData = "character",
				
				# === Set parameters for filtering and trimming ===			
				isPaired = "logical",
				filtMaxEE = "numeric",
				filtTruncQ = "numeric",
				filtTruncLen = "numeric",
				filtTrimLeft = "numeric",
				filtTrimRight = "numeric",
				filtMultiThread = "logical",
				filtVerbose = "logical",
				filtMatchIDs = "logical",
				filtMinLen = "numeric"
			),
			prototype = list(
				pathToData = "NA_character_",
				
				isPaired = FALSE,
				filtMaxEE = 2.5,
				filtTruncQ = 0, 
				filtTruncLen = 0,
				filtTrimLeft = 0,
				filtTrimRight = 0,
				filtMultiThread = TRUE, 
				filtVerbose = TRUE,
				filtMatchIDs = FALSE, 
				filtMinLen = 50
			),
			contains = "fastq2otu")

# Can extend and create with new("fastq2otu_single", ...)
#' Inherits methods and attributes from fastq2otu class
#' Lists paramaters that are essential for analysing paired-end data
#' @export
setClass("fastq2otu_paired",
			slots = c(
				# === Required inputs for paired-end data ===
				pathToData = "character",
				
				# === Set parameters for filtering and trimming ===			
				isPaired = "logical",
				filtMaxEE = "numeric",
				filtTruncQ = "numeric",
				filtTruncLen = "numeric",
				filtTrimLeft = "numeric",
				filtTrimRight = "numeric",
				filtMultiThread = "logical",
				filtVerbose = "logical",
				filtMatchIDs = "logical",
				filtMinLen = "numeric",

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
				
				isPaired = TRUE,
				filtMaxEE = c(2.5, 2.5),
				filtTruncQ = c(0, 0), 
				filtTruncLen = c(0, 0),
				filtTrimLeft = c(0, 0),
				filtTrimRight = c(0, 0),
				filtMultiThread = TRUE, 
				filtVerbose = TRUE,
				filtMatchIDs = TRUE, 
				filtMinLen = c(50, 50),
							
				doMergeSeqPairs = FALSE,
				mergeSeqPairsTrimOverhang = FALSE,
				mergeSeqPairsMinOverlap = 12,
				mergeSeqPairsMaxMismatch = 0,
				mergeSeqPairsReturnRejects = FALSE,
				mergeSeqPairsJustConcatenate = FALSE,
				mergeSeqPairsVerbose = FALSE
			),
			contains = "fastq2otu")

#' Directly create a fastq2otu object from paired-end data
#' 
#' Not all parameters can be modified using this method (i.e. verbose and multithread options are TRUE by default). 
#' @param database Required input. Path to fasta-formatted database containing reference sequences used for taxonomic assignment. 
#' @param in_dir Required input. Path to directory containing input FASTQ files. Forward and reverse reads must be in same directory.
#' @param out_dir Required input. Path to directory that will store all output files and directories. This directory may or may not exist.
#'
#' @param prefix Default is "myproject". Project-specific label that will be included in the names of all files and directories created.
#' @param id_list Default is NA. Path to new-line delimited text file containing SRA ids. This file will be used to run fastq-dump. 
#' @param run_fastq_dump Provide TRUE or FALSE. Default is FALSE. Determines whether Fastq-dump will be ran. id.list, get_fastq_script and fastq_dump 
#' parameters cannot be left as NA if TRUE. If FALSE, inputs for id.list, get_fastq_script and fastq_dump are ignored. 
#' @param fastq_dump Default is NA. Path to executable fastq-dump file. 
#' @param get_fastq_script Default is NA. Path to bash script `retrieve_sra_data.sh`, which allows fastq-dump script to be executed based on SRA ids provided in 
#' input file.
#' 
#' @param trim_primers Provide TRUE or FALSE. Default is FALSE. Determines whether bbduk.sh script will be ran. adapter.list, raw_data and 
#' bbduk_script parameters cannot be left empty if TRUE. DADA2 trimRight and trimLeft parameters can also be used to remove adapters of uniform 
#' length. The in.dir parameter should specify the path to output directory for the trimmed sequences while the raw_data parameter should specify 
#' the path to the untrimmed sequences. 
#' @param adapter_list Path to new-line delimited text file containing adapter sequences. This file will be used to run bbduk.sh.
#' @param raw_data Path to directory containing FASTQ-files with primers.
#' @param bbduk_script Path to executable bbduk.sh script. 
#'
#' @param max_err Sets maxEE parameter for DADA2 filterAndTrim function. Two numeric inputs are expected. 
#' @param trunc_qual Sets truncQ parameter for DADA2 filterAndTrim function. Two numeric inputs are expected.
#' @param trim_length Sets truncLen parameter for DADA2 filterAndTrim function. Two numeric inputs are expected.
#' @param trim_left Sets trimLeft parameter for DADA2 filterAndTrim function. Two numeric inputs are expected.
#' @param trim_right Sets trimRight parameter for DADA2 filterAndTrim function. Two numeric inputs are expected.
#' @param match_ids Sets matchIDs parameter for DADA2 filterAndTrim function. Two numeric inputs are expected.
#' @param min_length Sets minLen parameter for DADA2 filterAndTrim function. Two numeric inputs are expected.
#'
#' @param merge_pairs Provide TRUE or FALSE. Default is FALSE. Enables DADA2's mergePairs function to be executed. 
#' @param min_overlap Sets minOverlap parameter for DADA2 mergePairs function. 
#' @param max_mismatch Sets maxMismatch parameter for DADA2 mergePairs function. 
#' @param return_rejects Sets returnRejects parameter for DADA2 mergePairs function. 
#' @param just_concatenate Sets justConcatenate parameter for DADA2 mergePairs function. 
#' 
#' @param run_fastqc Provide TRUE or FALSE. Default is FALSE (parameter for fastqc_out will be ignored). Determines whether Fastqc will be executed. If TRUE, a value for fastqc_out must be provided.
#' @param fastqc_out Path to output directory that FASTQC can write to. 
#'
#' @param final_summary_table c 
#' @param plot_quality Provide TRUE or FALSE for creating aggregate DADA2 quality plots. 
#' @param plot_pdf Path to output PDF file containing DADA2 quality plots. 
#' @param err_pdf Path to output PDF file containing DADA2 learned error plots. 
#' @param dada_bandsize Sets BAND_SIZE parameter for DADA2's dada function.
#' @param dada_omegaA Sets OMEGA_A parameter for DADA2's dada function. 
#' @param chimera_detection_table Provide TRUE or FALSE for creating table showing all sequences labled as chimeric by DADA2 (uses DADA2's isBimeraDevovoTable function as opposed to removeBimeraDenovo). 
#' @param chimera_min_detection Sets minDetection parameter for DADA2's removeBimeraDenovo and isBimera functions. 
#' @param chimera_ignore_negatives Sets ignoreNegatives parameter for DADA2's removeBimeraDenovo and isBimera functions. 
#' @param chimera_minfold_parent_overabundance Sets minfoldParentOverabundance parameter for DADA2's removeBimeraDenovo and isBimera functions. 
#' @param chimera_dectection_parent_abundance Sets detectionParentAbundace parameter for DADA2's removeBimeraDenovo and isBimera functions. 			 
#' @param chimera_allow_one_off	Sets allowOneOff parameter for DADA2's removeBimeraDenovo and isBimera functions. 	
#' @param chimera_max_shift Sets maxShift parameter for DADA2's removeBimeraDenovo and isBimera functions. 
#' 
#' @param assign_tax_min_bootstrap Sets min-bootstrap value for DADA2 assignTaxonomy function.
#' @param assign_tax_try_complement	Provide TRUE or FALSE for tryComplement parameter for DADA2's assignTaxonomy function. 
#' @param assign_tax_levels Provide character list of taxonomic levels to include in taxonomy assignment. The default OTU table generated from this function includes assignments from kingdom to species.
#' 
#' @param merge_samples	Provide TRUE or FALSE to generate a merged sample table. Allows for cross-sample analyses, however depending on number of samples provided this tables may become very large.
#' @param final_merged_table Provide the name of the CSV file used to save the final merged table. If merge_samples is FALSE, this parameter can be ignored. 
#  @return S4 object of type fastq2otu_paired
#' @export
#' 
fastq2otu_paired <- function(in_dir, database, out_dir, isPaired = TRUE, prefix = "myproject",
							run_fastq_dump = FALSE, id_list = NA, trim_primers = FALSE, adapter.list = NA, raw_data = NA, 
							max_err = c(2.5, 2.5), trunc_qual = c(0, 0), trim_length = c(0, 0), trim_left = c(0, 0), trim_right = c(0, 0), min_length = c(50, 50), 
							merge_pairs = FALSE, min_overlap = 12 , max_mismatch = 0, return_rejects = FALSE, just_concatenate = FALSE, 
							run_fastqc = NA, fastqc_out = NA, final_summary_table = "seq_filter_log.txt", plot_quality = TRUE, plot_pdf = "quality_plots.pdf", err_pdf = "learn_errors_plots.pdf", dada_bandsize =16,
							dada_omegaA = 1e-10, chimera_dectection_table = FALSE, chimera_min_detection = 0.9, chimera_ignore_negatives = 1, chimera_minfold_parent_overabundance = 1.5, 
							chimera_dectection_parent_abundance = 2, chimera_allow_one_off = FALSE, chimera_max_shift = 4, 
							assign_tax_levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), assign_tax_min_bootstrap = 50, assign_tax_try_complement = TRUE,
							merge_samples = TRUE, final_merged_table = "final_merged_table.csv") {
							
	temp <- methods::new("fastq2otu_paired", 
		taxDatabase = database,
		pathToSingleData = in_dir,
		projectPrefix = prefix,
		outDir = out_dir,

		pathToSampleIDs = id_list,
		runFastqDump = run_fastq_dump,

		trimPrimers = trim_primers,
		listOfAdapters = adapter_list,
		pathToRawFastq = raw_data,

		filtMaxEE = max_err,
		filtTruncQ = trunc_qual,
		filtTruncLen = trim_length,
		filtTrimLeft = trim_left,
		filtTrimRight = trim_right,
		filtMinLen = min_length,
		
		doMergeSeqPairs = merge_pairs,
		mergeSeqPairsMinOverlap = min_overlap,
		mergeSeqPairsMaxMismatch = max_mismatch,
		mergeSeqPairsReturnRejects = return_rejects,
		mergeSeqPairsJustConcatenate = just_concatenate,
		
		runFastqc = run_fastqc
		pathToFastqcResults = fastqc_out,
	
		finalSummaryTable = final_summary_table,

		plotQuality = plot_quality,
		qualityPlotPDF = plot_pdf,
		learnErrorsPlotPDF = err_pdf,

		dadaBandSize = dada_bandsize,
		dadaOmegaA = dada_omegaA,

		createChimeraDetectionTable = chimera_dectection_table,
		chimeraDetectionMinSampleFraction = chimera_min_detection,
		chimeraDetectionIgnoreNegatives = chimera_ignore_negatives,
		chimeraDetectionMinFoldParentOverabundance = chimera_minfold_parent_overabundance,
		chimeraDetectionParentAbundance = chimera_dectection_parent_abundance,
		chimeraDetectionAllowOneOff = chimera_allow_one_off,
		chimeraDetectionMaxShift = chimera_max_shift,

		assignTaxMinBootstrap = assign_tax_min_bootstrap,
		assignTaxTryComplement = assign_tax_try_complement,
		assignTaxLevels = assign_tax_levels,

		mergeSamples = merge_samples,
		finalMergedTable = final_merged_table)
		
	return(temp)
}



#' Directly create a fastq2otu object from single-end data
#' 
#' Not all parameters can be modified using this method (i.e. verbose and multithread options are TRUE by default). 
#' @param database Required input. Path to fasta-formatted database containing reference sequences used for taxonomic assignment. 
#' @param in_dir Required input. Path to directory containing input FASTQ files. Forward and reverse reads must be in same directory.
#' @param out_dir Required input. Path to directory that will store all output files and directories. This directory may or may not exist.
#'
#' @param prefix Default is "myproject". Project-specific label that will be included in the names of all files and directories created.
#' @param id_list Default is NA. Path to new-line delimited text file containing SRA ids. This file will be used to run fastq-dump. 
#' @param run_fastq_dump Provide TRUE or FALSE. Default is FALSE. Determines whether Fastq-dump will be ran. id.list, get_fastq_script and fastq_dump 
#' parameters cannot be left as NA if TRUE. If FALSE, inputs for id.list, get_fastq_script and fastq_dump are ignored. 
#' @param fastq_dump Default is NA. Path to executable fastq-dump file. 
#' @param get_fastq_script Default is NA. Path to bash script `retrieve_sra_data.sh`, which allows fastq-dump script to be executed based on SRA ids provided in 
#' input file.
#' 
#' @param trim_primers Provide TRUE or FALSE. Default is FALSE. Determines whether bbduk.sh script will be ran. adapter.list, raw_data and 
#' bbduk_script parameters cannot be left empty if TRUE. DADA2 trimRight and trimLeft parameters can also be used to remove adapters of uniform 
#' length. The in.dir parameter should specify the path to output directory for the trimmed sequences while the raw_data parameter should specify 
#' the path to the untrimmed sequences. 
#' @param adapter_list Path to new-line delimited text file containing adapter sequences. This file will be used to run bbduk.sh.
#' @param raw_data Path to directory containing FASTQ-files with primers.
#' @param bbduk_script Path to executable bbduk.sh script. 
#'
#' @param max_err Sets maxEE parameter for DADA2 filterAndTrim function. One numeric input is expected. 
#' @param trunc_qual Sets truncQ parameter for DADA2 filterAndTrim function. One numeric input is expected. 
#' @param trim_length Sets truncLen parameter for DADA2 filterAndTrim function. One numeric input is expected. 
#' @param trim_left Sets trimLeft parameter for DADA2 filterAndTrim function. One numeric input is expected. 
#' @param trim_right Sets trimRight parameter for DADA2 filterAndTrim function. One numeric input is expected. 
#' @param match_ids Sets matchIDs parameter for DADA2 filterAndTrim function. One numeric input is expected. 
#' @param min_length Sets minLen parameter for DADA2 filterAndTrim function. One numeric input is expected. 
#'
#' 
#' @param run_fastqc Provide TRUE or FALSE. Default is FALSE (parameter for fastqc_out will be ignored). Determines whether Fastqc will be executed. If TRUE, a value for fastqc_out must be provided.
#' @param fastqc_out Path to output directory that FASTQC can write to. 
#'
#' @param final_summary_table c 
#' @param plot_quality Provide TRUE or FALSE for creating aggregate DADA2 quality plots. 
#' @param plot_pdf Path to output PDF file containing DADA2 quality plots. 
#' @param err_pdf Path to output PDF file containing DADA2 learned error plots. 
#' @param dada_bandsize Sets BAND_SIZE parameter for DADA2's dada function.
#' @param dada_omegaA Sets OMEGA_A parameter for DADA2's dada function. 
#' @param chimera_detection_table Provide TRUE or FALSE for creating table showing all sequences labled as chimeric by DADA2 (uses DADA2's isBimeraDevovoTable function as opposed to removeBimeraDenovo). 
#' @param chimera_min_detection Sets minDetection parameter for DADA2's removeBimeraDenovo and isBimera functions. 
#' @param chimera_ignore_negatives Sets ignoreNegatives parameter for DADA2's removeBimeraDenovo and isBimera functions. 
#' @param chimera_minfold_parent_overabundance Sets minfoldParentOverabundance parameter for DADA2's removeBimeraDenovo and isBimera functions. 
#' @param chimera_dectection_parent_abundance Sets detectionParentAbundace parameter for DADA2's removeBimeraDenovo and isBimera functions. 			 
#' @param chimera_allow_one_off	Sets allowOneOff parameter for DADA2's removeBimeraDenovo and isBimera functions. 	
#' @param chimera_max_shift Sets maxShift parameter for DADA2's removeBimeraDenovo and isBimera functions. 
#' 
#' @param assign_tax_min_bootstrap Sets min-bootstrap value for DADA2 assignTaxonomy function.
#' @param assign_tax_try_complement	Provide TRUE or FALSE for tryComplement parameter for DADA2's assignTaxonomy function. 
#' @param assign_tax_levels Provide character list of taxonomic levels to include in taxonomy assignment. The default OTU table generated from this function includes assignments from kingdom to species.
#' 
#' @param merge_samples	Provide TRUE or FALSE to generate a merged sample table. Allows for cross-sample analyses, however depending on number of samples provided this tables may become very large.
#' @param final_merged_table Provide the name of the CSV file used to save the final merged table. If merge_samples is FALSE, this parameter can be ignored. 
#'
#  @return S4 object of type fastq2otu_single
#' @export
fastq2otu_single <- function(function(in_dir, database, out_dir, isPaired = TRUE, prefix = "myproject",
							run_fastq_dump = FALSE, id_list = NA, trim_primers = FALSE, adapter_list = NA, raw_data = NA, 
							max_err = 2.5, trunc_qual = 0, trim_length = 0, trim_left = 0, trim_right = 0, min_length = 50, 
							run_fastqc = NA, fastqc_out = NA, final_summary_table = "seq_filter_log.txt", plot_quality = TRUE, plot_pdf = "quality_plots.pdf", err_pdf = "learn_errors_plots.pdf", 
							dada_bandsize =16, dada_omegaA = 1e-10, chimera_dectection_table = FALSE, chimera_min_detection = 0.9, chimera_ignore_negatives = 1, chimera_minfold_parent_overabundance = 1.5, 
							chimera_dectection_parent_abundance = 2, chimera_allow_one_off = FALSE, chimera_max_shift = 4, 
							assign_tax_levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), assign_tax_min_bootstrap = 50, assign_tax_try_complement = TRUE,
							merge_samples = TRUE, final_merged_table = "final_merged_table.csv") {
							
	temp <- methods::new("fastq2otu_single", 
		taxDatabase = database,
		pathToSingleData = in_dir,
		projectPrefix = prefix,
		outDir = out_dir,

		pathToSampleIDs = id_list,
		runFastqDump = fastq_dump,
		retrieveSRAData = get_fastq_script,

		trimPrimers = trim_primers,
		listOfAdapters = adapter_list,
		pathToRawFastq = raw_data,
		pathToUseBBDuk = bbduk_script,

		filtMaxEE = max_err,
		filtTruncQ = trunc_qual,
		filtTruncLen = trim_length,
		filtTrimLeft = trim_left,
		filtTrimRight = trim_right,
		filtMinLen = min_length,

		runFastqc = run_fastqc
		pathToFastqcResults = fastqc_out,
	
		finalSummaryTable = final_summary_table,

		plotQuality = plot_quality,
		qualityPlotPDF = plot_pdf,
		learnErrorsPlotPDF = err_pdf,

		dadaBandSize = dada_bandsize,
		dadaOmegaA = dada_omegaA,

		createChimeraDetectionTable = chimera_dectection_table,
		chimeraDetectionMinSampleFraction = chimera_min_detection,
		chimeraDetectionIgnoreNegatives = chimera_ignore_negatives,
		chimeraDetectionMinFoldParentOverabundance = chimera_minfold_parent_overabundance,
		chimeraDetectionParentAbundance = chimera_dectection_parent_abundance,
		chimeraDetectionAllowOneOff = chimera_allow_one_off,
		chimeraDetectionMaxShift = chimera_max_shift,

		assignTaxMinBootstrap = assign_tax_min_bootstrap,
		assignTaxTryComplement = assign_tax_try_complement,
		assignTaxLevels = assign_tax_levels,

		mergeSamples = merge_samples,
		finalMergedTable = final_merged_table)
		
	return(temp)
	
}

