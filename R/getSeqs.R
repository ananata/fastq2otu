#' Download SRA data using fastq-dump
#' Parses list of SRA accession numbers and uses fastq-dump to download FASTQ files directly from NCBI.
#' Downloaded files are gziped by default. 
#' 
#' @param object S4 object of type fastq2otu (i.e. fastq2otu_paired or fastq2otu_single)
#' @return path to directory
#' @export
getSeqs <- function(object) {
  # Extract required inputs
  fDumpScript <- system.file("exec", "fastq_dump", package = "fastq2otu")
  retrieveSRAData <- system.file("exec", "retrieve_sra_sequences.sh", package = "fastq2otu")
  sraList <- object@pathToSampleIDs
  output <- object@pathToData # Writes to input directory

  # Check to see if path exists
  if (!dir.exists(output)) {
    dir.create(output)
  }

  # Set system command
  command <- paste(retrieveSRAData, fDumpScript, sraList, output, sep = " ")

  # Execute command
  system(command)

  # Return string
  return(output)
}

