plotQuality <- function(configFile) {
  # Read YML file
  options <- yaml.load_file(configFile)

  # Extract required parameters
  fp <- options$pathToData
  out <- options$outDir
  label <- options$qualityPlotPDF

  # Verify that inputs are valid
  if (!dir.exists(fp) | length(list.files(path = fp, pattern = ".fastq")) == 0) {
    stop(sprintf("'%s' does not exist or no FASTQ files could be detected in path", fp))
  }

  if (!dir.exists(out)) {
    stop(sprintf("'%s' does not exist.", out))
  }

  # Extract all FASTQ files from path
  Fs <- sort(list.files(fp, pattern="*.fastq", full.names = TRUE))

  # Find current date
  currDate <- gsub("-", "", Sys.Date())

  ## Plot aggregate graph
  plotAgg <- dada2::plotQualityProfile(Fs, n = 1e+05, aggregate = T)

  ## Change plot title
  plotAgg <- plotAgg + ggplot2::ggtitle(paste(label, " Aggregate Quality Plot"))

  ## Return plot (only one plot is generated per dataset)
  return(plotAgg)

}

