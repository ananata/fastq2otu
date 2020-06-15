runFastqc <- function(path) {
  # Load required libraries
  require(yaml)
  require(fastqcr)
  require(dplyr)

  # Extract necessary parameters
  output <- options$pathToFastqcResults
  nThreads <- options$fastqcThreads
  experiment <- options$fastqcExperimentDescription
  label <- options$projectPrefix

  # Run FASTQC to obtain results
  fastqcr::fastqc(fq.dir = path, qc.dir = output, threads = nThreads)

  # Create report
  fastqcr::qc_report(output, result.file = file.path(output, paste0(label, "_fastqc_report")), experiment = experiment, preview = FALSE)

  # Return path to FASTQC results
  return(output)
}

