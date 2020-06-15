plotSeqErrors <- function(errObj) {
  # Read YML file (get configFile object from main.R)
  options <- yaml.load_file(configFile)

  # Extract required parameters
  doPlot <- options$plotErrors
  pdfF <- options$errPDF

  # Execute the function
  errPlot <- dada2::plotErrors(errObj)

  # Return plot
  return(errPlot)
}

