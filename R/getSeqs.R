#' Download SRA data using fastq-dump
#' Parses list of SRA accession numbers and uses fastq-dump to download FASTQ files directly from NCBI.
#' Downloaded files are gziped by default. 
#' 
#' @param object S4 object of type fastq2otu (i.e. fastq2otu_paired or fastq2otu_single)
#' @param useFastqDump Default is FALSE. If TRUE, fastqdump is used to download sequence files. 
#' @return path to directory
#' @export
getSeqs <- function(object, useFastqDump = FALSE) {
  # Get bash script
  # retrieveSRAData <- system.file("bash", "retrieve_sra_sequences.sh", package = "fastq2otu") -- Use after package is installed
  retrieveSRAData <- "/gpfs_fs/vahmp/users/Fettweis/Atopobium_NTA/Wright_Collaboration/r_package/fastq2otu/inst/bash/retrieve_sra_sequences.sh"

  if (useFastqDump) {
    # Extract required inputs
    # Get path to fastq-dump executable script
    # Check for objects returning NA/NULL or empty strings
    first_check <- ifelse(is.na(object@pathToFastqDump) | object@pathToFastqDump == "", TRUE, FALSE)
    if (first_check) { stop("Invalid path supplied for pathToFastqDump. Please check config file") }

    if (!is.null(object@pathToFastqDump) & file.exists(object@pathToFastqDump)) {
      fDumpScript <- object@pathToFastqDump
    } else {
      stop(paste0("Invalid input for pathToFastqDump. ", object@pathToFastqDump, " could not be found."))
    }

    # Get path to output directory
    # Check for objects returning NA/NULL or empty strings
    second_check <- ifelse(is.na(object@outDir) | object@outDir == "", TRUE, FALSE)
    if (second_check) { stop("Invalid path supplied for pathToData: ", object@outDir, ". Please check config file") }
    
    # Get path to sample IDs
    # Check for objects returning NA/NULL or empty strings
    ifelse(is.na(object@pathToSampleIDs) | object@pathToSampleIDs == "", stop("Invalid path supplied for pathToSampleIDs. Please check config file"), FALSE)
    if (!is.null(object@pathToSampleIDs) & file.exists(object@pathToSampleIDs)) {
      sraList <- object@pathToSampleIDs
      output <- object@outDir # Writes to input directory (validated when object was created)
    } else {
      stop(paste0("Invalid input for pathToSampleIDs. ", object@pathToSampleIDs, " could not be found."))
    }
    
    # Set system command - Example: ./retrieve_sra_sequences.sh ~/.sra-toolkit/bin/fastq-dump ../inst/examples/paired/paired-example_SRR_Acc_List.txt ./
    command <- paste(retrieveSRAData, fDumpScript, sraList, output, sep = " ")
    message("Executed: ", command, "\n")

  } else {
    # Get path to output directory
    # Check for objects returning NA/NULL or empty strings
    second_check <- ifelse(is.na(object@outDir) | object@outDir == "", TRUE, FALSE)
    if (second_check) { stop("Invalid path supplied for pathToData: ", object@outDir, ". Please check config file") }
     
    # Create directory if it does not exist
    if (!dir.exists(object@outDir)) {
	response <- readline(prompt=paste0("Are you sure you want to create ", object@outDir, "? <y/N> : "))
	if (response %in% c("Yes", "Y", "y", "yes")) {
 	  dir.create(object@outDir)
	  message("Created ", object@outDir)
	} else {
	  stop("Program was stopped. Could not create output directory")
	}
    }
    
    # Extract required inputs
    # Get path to sample URLs -- obtained from SRA Explorer site
    # Check for objects returning NA/NULL or empty strings
    ifelse(is.na(object@pathToSampleURLs) | object@pathToSampleURLs == "", stop("Invalid path supplied for pathToSampleURLs. Please check config file"), FALSE)
    if (!is.null(object@pathToSampleIDs) & file.exists(object@pathToSampleURLs)) {
      sample.urls <- object@pathToSampleURLs
    } else {
      stop(paste0("Invalid input for pathToSampleURLs. ", object@pathToSampleURLs, " could not be found."))
    }
    
    # Set system command
    if(.Platform$OS.type == "unix") {
      command <- paste0("wget -i ", sample.urls, " -P ", object@outDir)
      message("This is my output path: ", object@outDir)
      output <- object@outDir # Writes to input directory 
    } else {
	# TODO: Test on windows machine
      command <- paste0("xargs -n 1 curl -O < ", sample.urls)
      message("Executed: ", command, "\n")
      output <- object@outDir # Writes to input directory
    }
  }
  
  # Execute command
  system(command)

  # Return vector containing output directory and input command
  return(c(output, command))
}

root <- "/gpfs_fs/vahmp/users/Fettweis/Atopobium_NTA/Wright_Collaboration/r_package/fastq2otu/"

source(paste0(root, "R/readConfig.R"))
source(paste0(root, "R/setup.R"))

paired_config <- paste0(root, "inst/examples/paired/my_paired-example_config.yml")
paired_options <- yaml::yaml.load_file(paired_config)

object <- readConfig(paired_config, type = "seqdump")
fp <- getSeqs(object, useFastqDump = FALSE)

fp
