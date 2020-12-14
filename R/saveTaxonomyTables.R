#' Save OTU Tables
#' 
#' Saves list of tables as tab-delimited CSV files
#' @param tables List of  OTU tables
#' @param sample.names List of file ids
#' @return Path(s) to saved OTU tables
#' @export
saveTaxonomyTables <- function(table, sample.name, output.dir, index, label) {
  if (!dir.exists(output.dir)) {
	stop("Invalid path provided for output.dir")
  }

  # Create a directory for OTU Tables if they do not already exist
  if (!dir.exists(file.path(output.dir, paste0(label, "_taxonomy_tables")))) {
    dir.create(file.path(output.dir, paste0(label, "_taxonomy_tables")))
  }

   # Reorder columns
   Sequences <- rownames(table)
   new.table <- cbind(table, Sequences)
   rownames(new.table) <- NULL

   # Set .csv file name
   f.print <- file.path(output.dir, paste0(label, "_taxonomy_tables"), paste0(index, "_", sample.name, "_OTU_Table.csv"))
   rds.print <- file.path(output.dir, paste0(label, "_taxonomy_tables"), paste0(index, "_", sample.name, "_OTU_Table.rds"))

   # Save file
   write.table(new.table, file = f.print, sep = "\t")
   saveRDS(new.table, file = rds.print)

   # Return paths
   return(rds.print)
}

