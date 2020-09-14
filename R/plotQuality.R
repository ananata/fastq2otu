#" Generate Quality Distribution as Aggregrated Plot
#' 
#' @param fp Path to directory containing FASTQ files
#' @param out Ouput directory
#' @param label Unique file or project id used to generate file/plot title 
#' @export
plotQuality <- function(fp, out, label, object) {
	# Verify that inputs are valid
	if (!dir.exists(fp) | length(list.files(path = fp, pattern = ".fastq")) == 0) {
		stop(sprintf("'%s' does not exist or no FASTQ files could be detected in path", fp))
	}

	if (!dir.exists(out)) {
		stop(sprintf("'%s' does not exist.", out))
	}

	# Extract all FASTQ files from path
	Fs <- sort(list.files(fp, pattern="*.fastq", full.names = TRUE))

	# Find current date
	currDate <- gsub("-", "", Sys.Date())

	## Plot aggregate graph
	if (object@aggregate) {
	  plotAgg <- dada2::plotQualityProfile(Fs, n = object@qualN, aggregate = T)
	} else {
	  plotAgg <- lapply(Fs, plot_quality, n = object@qualN) # Returns a list of ggplots
	}
	if (!missing(label) & object@aggregate) {
		## Change plot title
		plotAgg <- plotAgg + ggplot2::ggtitle(paste0(label, " Aggregate Quality Plot"))
	} else if (!missing(label) & !object@aggregate) {
		plotAgg <- plotAgg + ggplot2::ggtitle(paste0(label, " Quality Plots")		
	}

	## Return plot (only one plot is generated per dataset)
	return(plotAgg)

}

#' Create a function the plots quality graphs on individual pages within a .pdf file
#' @param data Path(s) to input files
#' @param n Sampling number
#' @export
plot_quality <- function(data, n) {
	dada2::plotQualityProfile(data, n = n, aggregate = F)
}


