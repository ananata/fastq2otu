#' Generate FASTQC Report
#' @param object fastq2otu-type class object
#' @param path Path to directory containing input fastq files
#' @export
runFastqc <- function(object, path) {
	# Extract necessary parameters
	output <- object@pathToFastqcResults
	nThreads <- object@fastqcThreads
	experiment <- object@fastqcExperimentDescription
	label <- object@projectPrefix
	fastqc_path <- object@pathToFastqc
	install_fastqc <- object@installFastqc

	if (install_fastqc & is.na(fastqc_path)) {
		# Install fastqc at "~/bin/FastQC/fastqc"
		fastqcr::fastqc_install()
		fastqc_path <- "~/bin/FastQC/fastqc"
		message("FASTQC was successfully installed at: " fastqc_path)
	}
	if (install_fastqc & !is.na(fastqc_path)) {
		# Install fastqc at user-specfied destination
		# Error message is returned if platform is not Unix/Linux based
		if(!(.Platform$OS.type == "unix")) {
			stop("Unix system (MAC OSX or Linux) required.")
		} else { 
			fastqcr::fastqc_install(dest.dir = fastqc_path) 
		}
		message("FASTQC was successfully installed at: " fastqc_path)
	}
	# Run FASTQC to obtain results
	fastqcr::fastqc(fq.dir = path, qc.dir = output, threads = nThreads, fastqc.path = fastqc_path)

	# Create report
	fastqcr::qc_report(output, result.file = file.path(output, paste0(label, "_fastqc_report")), experiment = experiment, preview = FALSE)

	# Return path to FASTQC results
	return(output)
}

