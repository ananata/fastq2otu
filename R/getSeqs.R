#' Download SRA data using fastq-dump
#' Parses list of SRA accession numbers and uses fastq-dump to download FASTQ files directly from NCBI.
#' Downloaded files are gziped by default. 
#' 
#' @importFrom curl curl_download
#' @param object S4 object of type fastq2otu (i.e. fastq2otu_paired or fastq2otu_single)
#' @param useFastqDump Default is FALSE. If TRUE, fastqdump is used to download sequence files. 
#' @return path to directory
#' @export
getSeqs <- function(object, useFastqDump = FALSE) {
  # Get bash script
  retrieveSRAData <- system.file("bash", "retrieve_sra_sequences.sh", package = "FASTQ2OTU")

  if (useFastqDump) {
    # Extract required inputs
    # Get path to fastq-dump executable script
    # Check for objects returning NA/NULL or empty strings
    first_check <- ifelse(is.na(object@pathToFastqDump) | object@pathToFastqDump == "", TRUE, FALSE)
    if (first_check) { stop("Path Not Found Error: Invalid path supplied for pathToFastqDump. Please check config file") }

    if (!is.null(object@pathToFastqDump) & file.exists(object@pathToFastqDump)) {
      fDumpScript <- object@pathToFastqDump
    } else {
      stop(paste0("Invalid input for pathToFastqDump. ", object@pathToFastqDump, " could not be found."))
    }

    # Get path to output directory
    # Check for objects returning NA/NULL or empty strings
    second_check <- ifelse(is.na(object@outDir) | object@outDir == "", TRUE, FALSE)
    if (second_check) { stop("Path Not Found Error: Invalid path supplied for pathToData: ", object@outDir, ". Please check config file") }
    
    # Get path to sample IDs
    # Check for objects returning NA/NULL or empty strings
    ifelse(is.na(object@pathToSampleIDs) | object@pathToSampleIDs == "", stop("Invalid path supplied for pathToSampleIDs. Please check config file"), FALSE)
    if (!is.null(object@pathToSampleIDs) & file.exists(object@pathToSampleIDs)) {
      sraList <- object@pathToSampleIDs
      output <- object@outDir # Writes to input directory (validated when object was created)
    } else {
      stop(paste0("Path Not Found Error: Invalid input for pathToSampleIDs. ", object@pathToSampleIDs, " could not be found."))
    }
    
    # Set system command - Example: ./retrieve_sra_sequences.sh ~/.sra-toolkit/bin/fastq-dump ../inst/examples/paired/paired-example_SRR_Acc_List.txt ./
    command <- paste(retrieveSRAData, fDumpScript, sraList, output, sep = " ")
    message("Executed: ", command, "\n")

  } else {
    # Get path to output directory
    # Check for objects returning NA/NULL or empty strings
    second_check <- ifelse(is.na(object@outDir) | object@outDir == "", TRUE, FALSE)
    if (second_check) { stop("Invalid path supplied for pathToData: ", object@outDir, ". Please check config file") }
     
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
      command <- paste0("wget -i ", sample.urls, " -P ", " -c ", object@outDir)
      message("This is my output path: ", object@outDir)
      output <- object@outDir # Writes to input directory 
    } else {
      con = file(sample.urls, "r")
      write("\n", file = sample.urls, append = TRUE) # Append newline to end of file
      while ( TRUE ) {
         line = readLines(con, n = 1)
         if ( grepl('^ftp', line) ) {
	   # Download file only if it does not already exist
	   if (!file.exists(file.path(object@outDir, basename(line)))) {
           	curl::curl_download(line, file.path(object@outDir, basename(line)))
	   }
         } else {
	   break
	 }
      }
      close(con)
      output <- object@outDir # Writes to input directory
      return(output)
    }
  }
  
  # Execute command
  system(command)

  # Return vector containing output directory and input command
  return(output)
}


