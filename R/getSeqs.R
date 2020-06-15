getSeqs <- function() {
  # Load required library
  require(yaml)

  # Extract required inputs
  fDumpScript <- options$pathToFastqDump
  retrieveSRAData <- options$retrieveSRAData
  sraList <- options$pathToSamples
  output <- options$pathToPairedData

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

