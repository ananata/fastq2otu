saveSeqs <- function(tables, sample.names) {
  # Initialize a character vector
  pathsToSeqs <- vector("character", length = length(tables))

  # Create a directory for sequences if it does not already exist
  if (!dir.exists(file.path(options$outDir, "sequence_tables"))) {
    dir.create(file.path(options$outDir, "sequence_tables"))
  }

  # Save RDS files
  i <- 1
  for (table in tables) {
     if (typeof(table) == "logical") {
       # Set .rds file name
       seq.print <- file.path(options$outDir, "sequence_tables", paste0(i, "_", sample.names[i], "_seqtab_chimeras.rds"))
     } else {
       seq.print <- file.path(options$outDir, "sequence_tables", paste0(i, "_", sample.names[i], "_seqtab.rds"))
     }
     # Save path in list
     pathsToSeqs[i] <- seq.print

     # Save file
     dada2::saveRDS(table, file = seq.print)

     # Increment counter
     i <- i + 1
  }

  # Return true
  return(pathsToSeqs)
}

