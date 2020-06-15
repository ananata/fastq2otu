saveTaxonomyTables <- function(tables, sample.names) {
  # Initialize list
  pathsToOTUs <- vector("character", length = length(tables))

  # Create a directory for OTU Tables if they do not already exist
  if (!dir.exists(file.path(options$outDir, "taxonomy_tables"))) {
    dir.create(file.path(options$outDir, "taxonomy_tables"))
  }

  # Save OTU tables
  i <- 1
  for (table in tables) {
     # Reorder columns
     Sequences <- rownames(table)
     new.table <- cbind(table, Sequences)
     rownames(new.table) <- NULL

     # Set .csv file name
     f.print <- file.path(options$outDir, "taxonomy_tables", paste0(i, "_", sample.names[i], "_OTU_Table.csv"))

     # Add path to list
     pathsToOTUs[i] <- f.print

     # Save file
     write.table(new.table, file = f.print)

     # Increment counter
     i <- i + 1
  }

  # Return paths
  return(pathsToOTUs)

}

