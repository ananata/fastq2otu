#' Merge OTU tables from different samples
#' @param otutabs List of path(s) to OTU tables
#' @param seqtabs List of path(s) to frequency tables
#' @param final.print Name of file to save. 
#' @return Merged table that can be used to complete cross sample comparisons
#' @importFrom gtools mixedsort
#' @importFrom dada2 mergeSequenceTables
#' @export
mergeSamples <- function(otutabs, seqtabs, label, taxLevels) {
  # Verify that sequence tables correspond to OTU tables
  otuLabels <- na.omit(as.vector(sapply(strsplit(basename(otutabs), '_OTU_Table.rds'), '[', 1)))
  seqLabels <- na.omit(as.vector(sapply(strsplit(basename(seqtabs), '_seqtab.rds'), '[', 1)))

  if (gtools::mixedsort(otuLabels) != gtools::mixedsort(seqLabels)) {
    stop("IDs must match")
  }

  # Read sequence tables into R
  clean.seqtabs <- gtools::mixedsort(seqtabs[seqtabs != ""])
  seqtab.list <- lapply(clean.seqtabs, readRDS)

  # Merge seqtab.list (produces WIDE matrix)
  mergedSeqs <- dada2::mergeSequenceTables(tables = seqtab.list)
 
  # Transpose
  mergedSeqs <- t(mergedSeqs)
  Sequences <- rownames(mergedSeqs)
  mergedSeqs <- cbind(Sequences, mergedSeqs)
  rownames(mergedSeqs) <- NULL

  
  # Read OTU tables into R
  clean.otutabs <- gtools::mixedsort(otutabs[otutabs != ""])
  otutab.list <- lapply(clean.otutabs, readRDS)

  # Merge otutab.list (warnings are fixed in next step) by sequences and lowest tax level
  mergedOTU <- Reduce(function(x, y, ...) merge(x, y, all = TRUE, ...), otutab.list) # Fully merge all tables (i.e. # i.e. merge(...merge(merge(otutab1, otutab2), otutab3) ...))
  byCols <- colnames(mergedOTU)

  # Warnings: Duplicated column names - fixed here
  names(mergedOTU)[-c(1:length(taxLevels))] <- as.vector(sapply(strsplit(basename(otuLabels), "_"), '[', 2))

  # Reorder rows in mergedOTU to correspond to row order in mergedSeqs
  mergedSeqs <- as.data.frame(mergedSeqs)
  mergedOTU <- mergedOTU[match(mergedSeqs$Sequences, mergedOTU$Sequences), ]

  # Verify that OTU and sequence tables can be merged across
  # Randomly selects row
  i <- sample(1:nrow(mergedOTU), 1)
  j <- sample(1:nrow(mergedOTU), 1)

  # Verify that random sequence found in mergedOTU corresponds sequence found in mergedSeqs
  s1 <- mergedSeqs$Sequences[i] == mergedOTU$Sequences[i]
  s2 <- mergedSeqs$Sequences[j] == mergedOTU$Sequences[j]

  # Check that both tests (s1 and s2) yielded 'TRUE'
  if (!all(s1, s2)) {
    return(FALSE)
  }

  # Merge tables (replace columns in mergedOTU with columns in mergedSeqs)
  final.tab <- mergedOTU[ , 1:length(byCols)]
  final.tab <- cbind(final.tab, mergedSeqs[match(mergedOTU$Sequences, mergedSeqs$Sequences), 2:ncol(mergedSeqs)])

  # Write the table to file
  final.print <- paste0(label, "_final_merged_table.txt")
  write.table(final.tab, file = final.print, sep = "\t")

  # Return table
  return(final.print)
}

