#' Merge OTU tables from different samples
#' @param otutabs List of path(s) to OTU tables
#' @param seqtabs List of path(s) to frequency tables
#' @param final.print Name of file to save. 
#' @return Merged table that can be used to complete cross sample comparisons
#' @export
mergeSamples <- function(otutabs, seqtabs, final.print = "final_merged_table.csv") {
  # Load required libraries
  # Verify that sequence tables correspond to OTU tables
  otuLabels <- as.vector(sapply(strsplit(basename(otutabs), '_OTU_Table.csv'), '[', 1))
  seqLabels <- as.vector(sapply(strsplit(basename(seqtabs), '_seqtab.rds'), '[', 1))

  if (mixedsort(otuLabels) != mixedsort(seqLabels)) {
    stop("IDs must match")
  }

  # Read sequence tables into R
  seqtab.list <- lapply(gtools::mixedsort(seqtabs), readRDS)

  # Merge seqtab.list (produces WIDE matrix)
  mergedSeqs <- dada2::mergeSequenceTables(tables = seqtab.list)

  # Transpose
  mergedSeqs <- t(mergedSeqs)
  Sequences <- rownames(mergedSeqs)
  mergedSeqs <- cbind(Sequences, mergedSeqs)
  rownames(mergedSeqs) <- NULL

  # Read OTU tables into R
  otutab.list <- lapply(gtools::mixedsort(otutabs), read.table)

  # Merge otutab.list (warnings are fixed in next step)
  byCols <- c(options$assignTaxLevels, "Sequences")
  mergedOTU <- Reduce(function(x, y) (merge(x, y, by = byCols, all = TRUE)), otutab.list,
                accumulate = FALSE) # i.e. merge(...merge(merge(otutab1, otutab2), otutab3) ...)

  # Warnings: Duplicated column names - fixed here
  names(mergedOTU)[-c(1:length(byCols))] <- as.vector(sapply(strsplit(basename(otuLabels), "_"), '[', 2))

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
    return("You can not merge across")
  }

  # Merge tables (replace columns in mergedOTU with columns in mergedSeqs)
  final.tab <- mergedOTU[ , 1:length(byCols)]
  final.tab <- cbind(final.tab, mergedSeqs[match(mergedOTU$Sequences, mergedSeqs$Sequences), 2:ncol(mergedSeqs)])

  # Write the table to file
  write.table(final.tab, file = final.print, sep = "\t")

  # Return table
  return(final.tab)
}

