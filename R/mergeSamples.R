#' Merge OTU tables from different samples
#' @param otutabs List of path(s) to OTU tables
#' @param seqtabs List of path(s) to frequency tables
#' @param final.print Name of file to save.
#' @return Merged table that can be used to complete cross sample comparisons
#' @importFrom gtools mixedsort
#' @importFrom dada2 mergeSequenceTables
#' @importFrom reshape2 melt
#' @importFrom data.table setDT
#' @export
mergeSamples <- function(otutabs, seqtabs, label, taxLevels) {
  # Verify that sequence tables correspond to OTU tables
  otuLabels <- na.omit(as.vector(sapply(strsplit(basename(otutabs), '_OTU_Table.rds'), '[', 1)))
  seqLabels <- na.omit(as.vector(sapply(strsplit(basename(seqtabs), '_seqtab.rds'), '[', 1)))
  
  if (gtools::mixedsort(otuLabels) != mixedsort(seqLabels)) {
    stop("IDs must match")
  }
  
  # Read sequence tables into R
  otutab.list <- lapply(mixedsort(otu.paths[otu.paths != ""]), readRDS)
  seqtab.list <- lapply(mixedsort(seq.paths[seq.paths != ""]), readRDS)
  
  mergedSeqs <- data.table::setDT(as.data.frame(t(dada2::mergeSequenceTables(tables = seqtab.list))), keep.rownames = "Sequences")
  mergedOTU <- Reduce(function(x, y, ...) merge(x, y, all = TRUE, ...), otutab.list)
  
  # Reorder rows in mergedOTU to correspond to row order in mergedSeqs
  ordered.mergedOTU <- mergedOTU[match(mergedSeqs$Sequences, mergedOTU$Sequences), ]
  
  # Combine columns
  final.tab <- merge(ordered.mergedOTU, mergedSeqs)
  
  # Create proportional table
  metaCols <- as.integer(length(taxLevels)) + 1
  prop.tab <- cbind(final.tab[ , 1:metaCols] , prop.table(as.matrix(final.tab[ , -c(1:metaCols)]), margin = 2))
  
  # Find maximum percent abundancies for each sequence variant in the table
  # Find average (factors in a lot of zeros)
  meanPercent <- apply(prop.tab[ , -c(1:metaCols)], 1, FUN=mean)
  maxPercent <- apply(prop.tab[ , -c(1:metaCols)], 1, FUN=max)
  
  # Rearrage columns
  final.prop.tab <- cbind(prop.tab[ , 1:metaCols], meanPercent, maxPercent, prop.tab[ , -c(1:metaCols)])
  
  # Write the table to file
  final.print <- paste0(label, "_final_merged_frequency_table.txt")
  f.rdata <- paste0(label, "_final_merged_frequency_table.rds")
  write.table(final.tab, file = final.print, sep = "\t")
  saveRDS(final.tab, file = f.rdata)
  
  prop.print <- paste0(label, "_final_merged_percent_table.txt")
  p.rdata <- paste0(label, "_final_merged_percent_table.rds")
  write.table(final.prop.tab, file = prop.print, sep = "\t")
  saveRDS(final.prop.tab, file = p.rdata, compress = T)
  
  # Return table
  return(final.print)
}

