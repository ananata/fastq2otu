#" Generate Quality Distribution as Aggregrated Plot
#' @title plotQuality
#' @param fp Path to directory containing FASTQ files
#' @param label Unique file or project id used to generate file/plot title
#' @param object fastPlotQuality object 
#' @export
plotQuality <- function(fp, label, object) {
	FLAG <- FALSE
	# Verify that inputs are valid
        if (length(fp) == 1 & dir.exists(fp)) {
                # Extracts all FASTQ from path
                Fs <- sort(list.files(fp, pattern="*.fastq(.gz)?", full.names = TRUE))
		if (length(Fs) < 1) {
			stop("No FASTQ files could be detected in path")
		}
        }	
	else if (length(fp) > 1) {
		Fs <- fp
		FLAG <- TRUE
	}  else {
		stop(sprintf("'%s' does not exist", fp))
	}

	## Plot aggregate graph
	if (object@aggregateQual) {
	  plotAgg <- dada2::plotQualityProfile(Fs, n = object@qualN, aggregate = T) # Returns a single plot
	} else {
	  plotAgg <- lapply(Fs, plot_quality, n = object@qualN) # Returns a list of ggplots
	}

	label <- object@projectPrefix
	if (!missing(label) & object@aggregateQual) {
		## Change plot title
		plotAgg <- plotAgg + ggplot2::ggtitle(paste0(label, " Aggregate Quality Plot"))
	} else if (!missing(label) & !object@aggregateQual) {
		plotAgg <- plotAgg + ggplot2::ggtitle(paste0(label, " Quality Plots"))		
	}

	## Return plot (only one plot is generated per dataset)
	if (FLAG) {
		save(plotAgg, file = "quality_distribution_AFTER.RData")
	} else {
		save(plotAgg, file = "quality_distribution_BEFORE.RData")
	}
	return(plotAgg)

}

#' Create a function the plots quality graphs on individual pages within a .pdf file
#' @param data Path(s) to input files
#' @param n Sampling number
#' @export
plot_quality <- function(data, n) {
	dada2::plotQualityProfile(data, n = n, aggregate = F)
}


