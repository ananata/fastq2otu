#' Download SRA data using fastq-dump
#' Parses list of SRA accession numbers and uses fastq-dump to download FASTQ files directly from NCBI.
#' Downloaded files are gziped by default. 
#' 
#' @param object S4 object of type fastq2otu (i.e. fastq2otu_paired or fastq2otu_single)
#' @return path to directory
#' @export
getSeqs <- function(object) {
  # Get bash script
  retrieveSRAData <- system.file("exec", "retrieve_sra_sequences.sh", package = "fastq2otu")
  
  if (object@useDump) {
    # Extract required inputs
    # Get path to fastq-dump executable script
    if (!is.null(object@pathToFastqDump) & file.exists(object@pathToFastqDump)) {
      fDumpScript <- object@pathToFastqDump
    } else {
      stop(paste0("Invalid input for pathToFastqDump.", object@pathToFastqDump, " could not be found."))
    }
    
    # Get path to sample IDs
    if (!is.null(object@pathToSampleIDs) & !file.exists(object@pathToSampleIDs)) {
      sraList <- object@pathToSampleIDs
      output <- object@pathToData # Writes to input directory (validated when object was created)
    }
    
    # Set system command
    command <- paste(retrieveSRAData, fDumpScript, sraList, output, sep = " ")
  } else {
    # Extract required inputs
    # Get path to sample URLs -- obtained from SRA Explorer site
    if (file.exists(object@pathToSampleURLs)) {
      sample.urls <- object@pathToSampleURLs
    } else {
      stop(paste0("Invalid input for pathToSampleURLs", object@pathToSampleURLs, " could not be found."))
    }
    
    # Set system command
    if(.Platform$OS.type == "unix") {
      command <- paste0("wget -i ", sample.urls, " -P ", object@pathToData) 
    } else {
      command <- paste0("xargs -n 1 curl -O < ", sample.urls)
      output <- object@pathToData # Writes to input directory
    }
  }
  
  # Execute command
  system(command)

  # Return string
  return(output)
}

