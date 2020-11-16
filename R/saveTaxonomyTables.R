#' Save OTU Tables
#' 
#' Saves list of tables as tab-delimited CSV files
#' @param tables List of  OTU tables
#' @param sample.names List of file ids
#' @return Path(s) to saved OTU tables
#' @export
saveTaxonomyTables <- function(tables, sample.names, output.dir, label) {
  # Initialize list
  pathsToOTUs <- vector("character", length = length(tables))

  # Create a directory for OTU Tables if they do not already exist
  if (!dir.exists(file.path(output.dir, paste0(label, "_taxonomy_tables")))) {
    dir.create(file.path(output.dir, paste0(label, "_taxonomy_tables")))
  }

  # Save OTU tables
  for (i in 1:length(sample.names)) {
     # Reorder columns
     Sequences <- rownames(tables[i])
     new.table <- cbind(tables[i], Sequences)
     rownames(new.table) <- NULL

     # Set .csv file name
     f.print <- file.path(output.dir, paste0(label, "_taxonomy_tables"), paste0(i, "_", sample.names[i], "_OTU_Table.csv"))

     # Add path to list
     pathsToOTUs[i] <- f.print

     # Save file
     write.table(new.table, file = f.print, sep = "\t")
  }

  # Return paths
  return(pathsToOTUs)
}

