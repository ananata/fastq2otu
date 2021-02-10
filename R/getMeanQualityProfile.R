#' Get quality profile of a fastq file(s).
#' 
#' This function generates a table that reports the mean quality scores
#' observed at each sequence position for the input fastq file(s).
#' 
#' Code modified from DADA2's plotQualityProfile function
#' 
#' 
#' @param fl (Required). \code{character}.
#'  File path(s) to fastq or fastq.gz file(s).
#' 
#' @param n (Required). Default 500,000.
#'  The number of records to sample from the fastq file.
#' 
#' @param project (Required). Provide a character label that describes the data being handled.  
#' 
#' @param layout (Required). Default is Single. 
#'  Provide "Single" if handling single-end data; "Forward" if handling foward-reads, 
#'  and "Reverse" if handling reverse-reads. 
#'  
#' @return A table with the following column names: c("Position", "Mean", "Orientation", "Study")
#'  
#' @importFrom ShortRead qa
#' 
#'
#' @example getMeanQualityProfile(f1, project = "Study1")
#' @export
#' 
getMeanQualityProfile <- function(fl, n=500000, project, layout = "Single") {
  aggregate=TRUE 
  statdf <- data.frame(Cycle=integer(0), Mean=numeric(0), file=character(0))
  anndf <- data.frame(minScore=numeric(0), label=character(0), rclabel=character(0), rc=numeric(0), file=character(0))
  
  FIRST <- TRUE
  for(f in fl[!is.na(fl)]) {
    srqa <- ShortRead::qa(f, n=n)
    df <- srqa[["perCycle"]]$quality
    rc <- sum(srqa[["readCounts"]]$read) # Handle aggregate form from qa of a directory
    if (rc >= n) { 
      rclabel <- paste("Reads >= ", n)
    } else {
      rclabel <- paste("Reads: ", rc)
    }
    # Calculate summary statistics at each position
    means <- rowsum(df$Score*df$Count, df$Cycle)/rowsum(df$Count, df$Cycle)

    if(FIRST) {
      plotdf <- cbind(df, file=basename(f))
      FIRST <- FALSE
    } else { plotdf <- rbind(plotdf, cbind(df, file=basename(f))) }
    statdf <- rbind(statdf, data.frame(Cycle=as.integer(rownames(means)), Mean=means, file=basename(f)))
    anndf <- rbind(anndf, data.frame(minScore=min(df$Score), label=basename(f), rclabel=rclabel, rc=rc, file=basename(f)))
  }
  # Create plot
  plotdf.summary <- aggregate(Count ~ Cycle + Score, plotdf, sum)
  plotdf.summary$label <- paste(nrow(anndf), "files (aggregated)")
  means <- rowsum(as.numeric(plotdf.summary$Score*plotdf.summary$Count), as.numeric(plotdf.summary$Cycle)/rowsum(plotdf.summary$Count), plotdf.summary$Cycle)
  statdf.summary <- data.frame(Position=as.integer(rownames(means)), Mean=means, Orientation=rep(layout, length(means)),
                               Study=rep(project, length(means)))
  return(statdf.summary)
}

