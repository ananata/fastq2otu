#' Save sequence frequency tables in RDS file
#'
#' @param tables Sequence frequency table.
#' @param sample.names SRA/fastq file ids (extracted directly from file)
#' @param output.dir Output directory to write RDS file
#' @param add.table Generates an additional table storing all chimeric sequences flagged by DADA2's 
#' @export 
saveSeqs <- function(table, sample.name, index,  output.dir, add.table = FALSE, label) {

  # Check to make sure all inputs are valid
  if (is.null(table)) {
    stop("No tables were found.")
  }

  if (is.null(sample.name)) {
    stop("No sample names were found.")
  }

  if (is.na(output.dir) |!dir.exists(output.dir)) {
    stop("Output directory could not be found")
  }

  # Create a directory for sequences if it does not already exist
  if (missing(label)) { label <- Sys.Date() }

  if (!dir.exists(file.path(output.dir, paste0(label, "_sequence_tables")))) {
    dir.create(file.path(output.dir, paste0(label, "_sequence_tables")))
  }

  # Save RDS files
  if (add.table == TRUE) {
       # Set .rds file name
       seq.print <- file.path(output.dir, paste0(label, "_sequence_tables"), paste0(index, "_", sample.name, "_seqtab_chimeras.rds"))
  } else {
       seq.print <- file.path(output.dir, paste0(label, "_sequence_tables"), paste0(index, "_", sample.name, "_seqtab.rds"))
  }
  
  # Save path in list
  pathsToSeq <- seq.print

  # Save file
  if (!is.na(table)) {
    rownames(table) <- sample.name
    base::saveRDS(table, file = seq.print)
  }

  # Return true
  return(pathsToSeq)
}

