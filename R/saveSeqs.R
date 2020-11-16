#' Save sequence frequency tables in RDS file
#'
#' @param tables List of sequence frequency tables
#' @param sample.names List of SRA/fastq file ids (extracted directly from files)
#' @param output.dir Output directory to write RDS files
#' @param add.table Generates a table storing all chimeric sequences flagged by DADA2's 
#' @export 
saveSeqs <- function(tables, sample.names, output.dir, add.table = FALSE, label) {
  # Initialize a character vector
  pathsToSeqs <- vector("character", length = length(tables))

  # Check to make sure all inputs are valid
  if (is.null(tables) | length(tables) == 0) {
    stop("No tables were found.")
  }

  if (is.null(sample.names) | length(sample.names) == 0) {
    stop("No sample names were found.")
  }

  if (is.na(output.dir)) {
    stop("No output directory was provided")
  }

  # Create a directory for sequences if it does not already exist
  if (!dir.exists(file.path(output.dir, paste0(label, "_sequence_tables")))) {
    dir.create(file.path(output.dir, paste0(label, "_sequence_tables")))
  }

  # Save RDS files
  for (i in 1:length(sample.names)) {
     if (add.table == TRUE) {
       # Set .rds file name
       seq.print <- file.path(output.dir, paste0(label, "_sequence_tables"), paste0(i, "_", sample.names[i], "_seqtab_chimeras.rds"))
     } else {
       seq.print <- file.path(output.dir, paste0(label, "_sequence_tables"), paste0(i, "_", sample.names[i], "_seqtab.rds"))
     }
     # Save path in list
     pathsToSeqs[i] <- seq.print

     # Save file
     if (!is.na(tables[i])) 
     	base::saveRDS(tables[i], file = seq.print)
  }

  # Return true
  return(pathsToSeqs)
}

