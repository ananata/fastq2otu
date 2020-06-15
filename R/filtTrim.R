# Update to include paired-end data
filtTrim <- function(Fs) {
  require(dada2)
  require(yaml)

  # Extract required inputs
  label <- options$projectPrefix
  path <- options$outDir


  # Extract SRA run numbers from filenames (i.e. SRR829719)
  sample.names <- sapply(strsplit(basename(Fs), ".fastq"), `[`, 1)

  # Create output file names for filtered sequences
  filtFs <- file.path(path, paste0(label, "_filtered"), paste0(sample.names, "_filt_trimmed.fastq"))

  # Check to see if files exist
  if (all(file.exists(filtFs))) {
    return(filtFs)
  } else {
    notExist <- subset(filtFs, !file.exists(filtFs))
  }

  # Filter and Trim (only executes paths that do not already exist)
  filtTrim <- filterAndTrim(Fs, notExist, maxEE=options$filtMaxEE, trimLeft=options$filtTrimLeft,
                            trimRight=options$filtTrimRight, truncLen=options$filtTruncLen, multithread=options$filtMultiThread,
                            verbose=options$filtVerbose, minLen=options$filtMinLen, matchIDs=options$filtMatchIDs, truncQ=options$filtTruncQ)


  # Save important R objects (needed for summary table)
  save(filtTrim, file = "filtered_objects.RData")

  # Return paths
  return(filtFs)
}

