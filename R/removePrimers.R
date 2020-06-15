removePrimers <- function() {
  # Extract required inputs
  adapters <- options$listOfAdapters
  samples <- options$listOfSamples
  inDir <- options$pathToRawFastq
  outDir <- options$pathToNoPrimers
  bbdukScript <- options$pathToUseBBDuk

  # Set system command
  command <- paste(bbdukScript, inDir, outDir, samples, adapters, sep = " ")

  # Execute command
  system(command)

  # Return string
  return(outDir)
}

