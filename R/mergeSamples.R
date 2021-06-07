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
  transposed.mergedSeqs <- t(mergedSeqs)
  Sequences <- rownames(transposed.mergedSeqs)
  transposed.mergedSeqs <- cbind(Sequences, transposed.mergedSeqs)
  rownames(transposed.mergedSeqs) <- NULL  
  
  # Read OTU tables into R
  clean.otutabs <- gtools::mixedsort(otutabs[otutabs != ""])
  otutab.list <- lapply(clean.otutabs, read.csv)

  # Merge otutab.list (warnings are fixed in next step) by sequences and lowest tax level
  mergedOTU <- Reduce(function(x, y, ...) merge(x, y, all = TRUE, ...), otutab.list) # Fully merge all tables (i.e. # i.e. merge(...merge(merge(otutab1, otutab2), otutab3) ...))
  byCols <- colnames(mergedOTU)

  # Warnings: Duplicated column names - fixed here
  names(mergedOTU)[-c(1:length(taxLevels))] <- as.vector(sapply(strsplit(basename(otuLabels), "_"), '[', 2))

  # Reorder rows in mergedOTU to correspond to row order in mergedSeqs
  transposed.mergedSeqs <- as.data.frame(transposed.mergedSeqs)
  mergedOTU <- mergedOTU[match(transposed.mergedSeqs$Sequences, mergedOTU$Sequences), ]
  
  # Create proportional table
  melted.mergedSeqs <- reshape2::melt(data.table::setDT(as.data.frame(mergedSeqs), keep.rownames = TRUE), "rn")

  # Create new column
  total <- as.matrix(by(melted.mergedSeqs[ , -c(1:2), drop=FALSE], melted.mergedSeqs$rn, FUN=colSums))

  # Transpose
  transposed.total <- t(total)
 
  # Copy final.tab
  # Verify that OTU and sequence tables can be merged across
  # Randomly selects row
  i <- sample(1:nrow(mergedOTU), 1)
  j <- sample(1:nrow(mergedOTU), 1)

  # Verify that random sequence found in mergedOTU corresponds sequence found in mergedSeqs
  s1 <- transposed.mergedSeqs$Sequences[i] == mergedOTU$Sequences[i]
  s2 <- transposed.mergedSeqs$Sequences[j] == mergedOTU$Sequences[j]

  # Check that both tests (s1 and s2) yielded 'TRUE'
  if (!all(s1, s2)) {
    return(FALSE)
  }

  # Merge tables (replace columns in mergedOTU with columns in mergedSeqs)
  final.tab <- mergedOTU[ , 1:length(byCols)]
  final.tab <- cbind(final.tab, transposed.mergedSeqs[match(mergedOTU$Sequences, transposed.mergedSeqs$Sequences), 2:ncol(transposed.mergedSeqs), drop = FALSE])

  # Copy the final table
  prop.tab <- final.tab

  # Iterate through both tables and divide column values by their totals
  # Obtain the percent abundance of sequences within their respective datasets
  for (tName in colnames(transposed.total)) { # Iterate through totals
        for (cName in colnames(prop.tab)[-c(1:length(byCols))]) { # Iterate through counts
                # Check to see if column names are identical
                if (tName == cName) {
                        # Calculate percents
                        prop.tab[ , cName] <- (as.numeric(prop.tab[ , cName]) / as.numeric(transposed.total[1, tName])) * 100
                }
        }
  }

  # Add statistical information
  # Find minimum and maximum percent abundancies for each row in the table
  # Find average (take with grain of salt - there are a lot of zeros)
  prop.tab$Avg_Percent <- apply(prop.tab[ , -c(1:length(byCols))], 1, FUN=mean)
  prop.tab$Min_Percent <- apply(prop.tab[ , -c(1:length(byCols))], 1, FUN=min)
  prop.tab$Max_Percent <- apply(prop.tab[ , -c(1:length(byCols))], 1, FUN=max)
  
  # Rearrage columns
  prop.tab <- prop.tab[ , c(1:length(byCols), (ncol(prop.tab) - 1), ncol(prop.tab), (ncol(prop.tab) - 2), (length(byCols) + 1):(ncol(prop.tab) - 3))]

  # Write the table to file
  final.print <- paste0(label, "_final_merged_frequency_table.txt")
  f.rdata <- paste0(label, "_final_merged_frequency_table.rds")
  write.table(final.tab, file = final.print, sep = "\t")
  saveRDS(final.tab, file = f.rdata)
  message("Created: ")
  message(final.print)
  message(f.rdata)

  prop.print <- paste0(label, "_final_merged_percent_table.txt")
  p.rdata <- paste0(label, "_final_merged_percent_table.rds")
  write.table(prop.tab, file = prop.print, sep = "\t")
  saveRDS(prop.tab, file = p.rdata)
  message("Created: ")
  message(prop.print)
  message(p.rdata)

  # Return table
  return(final.print)
}

