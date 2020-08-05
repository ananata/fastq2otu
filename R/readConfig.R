#' Read in a config file (preferred method)
#'
#' This function parses config file to extract data needed to create a FASTQ2OTU class object. It assumes that
#' the file is in YML format and contains default key terms.
#' Execute getSlots("fastq2otu") to get all required arguments or refer to documenation for more details and/or
#' examples.
#'
#' Use when config file is supplied and all required values are provided.
#' @param configFile Path to config file (YML-formatted)
#' @param isPaired Default is FALSE. If TRUE, parameters for paired-end data will be used. 
#' @param type Default is 'auto'. String or vector containing the following value(s): 'auto', 'assignTax', 'report', 'seqdump', 'primertrim', 'filter', 'qualityplot'. Used to specify what type of object is being created.
#' If a vector is provided it must contain one of the following combination of strings c('auto', 'assignTax') and  c('auto', 'filter')
#' @return An S4 object
#' @keyword internal
#' @export
readConfig <- function(configFile, isPaired = FALSE, type = 'auto') {
	if (!file.exists(configFile)) {
			stop(sprintf("%s could not be found, please enter a valid path", configFile))
	}
	options <- yaml::yaml.load_file(configFile)

	# Create new fastq2otu class object based on sequencing method
	if ((type == "assignTax" | "assignTax" %in% type) & isPaired) {
		if (length(type) == 1) {
			temp <- setFastAssignTaxa(refDatabase = options$taxDatabase, 
					prefix = ifelse(!is.null(options$projectPrefix), options$projectPrefix, "myproject"),
					minBootstrap = ifelse(!is.null(options$assignTaxMinBootstrap), options$assignTaxMinBootstrap, 50), 
					tryComplement = ifelse(!is.null(options$assignTaxTryComplement), options$assignTaxTryComplement, FALSE), 
					showBootstraps = ifelse(!is.null(options$assignTaxOutputBootstraps), options$assignTaxOutputBootstraps, FALSE), 
					taxLevels = ifelse(!is.null(options$assignTaxLevels), options$assignTaxLevels, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")),
					verbose = ifelse(!is.null(options$verbose), options$verbose, FALSE), 
					multithread = ifelse(!is.null(options$multithread), options$multithread, FALSE))
		} 
		else if ('auto' %in% type & length(type) == 2) {
			temp <- setFastAssignTaxa(inDir = options$pathToData,
						outDir = options$outDir, 
						mergeSeqs = ifelse(!is.null(options$mergePairs), options$mergePairs, FALSE), 
						trimOverhang = ifelse(!is.null(options$mergePairsTrimOverhang), options$mergePairsTrimOverhang, FALSE),
						minOverlap = ifelse(!is.null(options$mergePairsMinOverlap), options$mergePairsMinOverlap, 12),
						maxMismatch = ifelse(!is.null(options$mergePairsMaxMismatch), options$mergePairsMaxMismatch, 0), 
						returnRejects = ifelse(!is.null(options$mergePairsReturnRejects), options$mergePairsReturnRejects, FALSE), 
						justConcatenate = ifelse(!is.null(options$mergePairsJustConcatenate), options$mergePairsJustConcatenate, FALSE),
						verbose = ifelse(!is.null(options$verbose), options$verbose, FALSE),
						prefix = ifelse(!is.null(options$projectPrefix), options$projectPrefix, "myproject"),
						isPaired = TRUE, 
						derepN = ifelse(!is.null(options$derepN), as.numeric(as.character(options$derepN)), 1e+06),
						getErrPDF = ifelse(!is.null(options$saveErrorsPlot), options$saveErrorsPlot, FALSE),
						errN = ifelse(!is.null(options$errN), as.numeric(as.character(options$errN)), 1e+08), 
						multithread = ifelse(!is.null(options$multithread), options$multithread, FALSE),
						dadaBandSize = ifelse(!is.null(options$dadaBandSize), as.numeric(as.character(options$dadaBandSize)), 16), 
						dadaOmegaA = ifelse(!is.null(options$dadaOmegaA), as.numeric(as.character(options$dadaOmegaA)), 1e-40), 
						getChimeraTable = ifelse(!is.null(options$createChimeraDetectionTable), options$createChimeraDetectionTable, FALSE), 
						minSampleFraction = ifelse(!is.null(options$chimeraDetectionMinSampleFraction), as.numeric(as.character(options$chimeraDetectionMinSampleFraction)), 0.9), 
						ignoreNegatives = ifelse(!is.null(options$chimeraDetectionIgnoreNegatives), as.numeric(as.character(options$chimeraDetectionIgnoreNegatives)), 1), 
						minFoldParentOverAbundance = ifelse(!is.null(options$chimeraDetectionMinFoldParentOverabundance), as.numeric(as.character(options$chimeraDetectionMinFoldParentOverabundance)), 1.5), 
						detectionAbundance = ifelse(!is.null(options$chimeraDetectionParentAbundance), as.numeric(as.character(options$chimeraDetectionParentAbundance)), 2), 
						allowOneOff = ifelse(!is.null(options$chimeraDetectionAllowOneOff), options$chimeraDetectionAllowOneOff, FALSE), 
						maxShift = ifelse(!is.null(options$chimeraDetectionMaxShift), as.numeric(as.character(options$chimeraDetectionMaxShift)), 16),
						
						minBootstrap = ifelse(!is.null(options$assignTaxMinBootstrap), as.numeric(as.character(options$assignTaxMinBootstrap)), 50), 
						tryComplement = ifelse(!is.null(options$assignTaxTryComplement), as.numeric(as.character(options$assignTaxTryComplement)), FALSE), 
						showBootstraps = ifelse(!is.null(options$assignTaxOutputBootstraps), options$assignTaxOutputBootstraps, FALSE), 
						taxLevels = ifelse(!is.null(options$assignTaxLevels), options$assignTaxLevels, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")))
		}
		else {
			stop("Invalid input provided for type parameter")
		}
	}
	else if ((type == "filter" | "filter" %in% type) & isPaired) {
		if (length(type) == 1) {
			temp <- setFastFilter(inDir = options$pathToData,
					outDir = options$outDir,
					filtVerbose = ifelse(!is.null(options$verbose), options$verbose, FALSE),
					prefix = ifelse(!is.null(options$projectPrefix), options$projectPrefix, "myproject"), 
					maxEE = ifelse(!is.null(options$filtMaxEE), as.numeric(as.character(options$filtMaxEE)), c(2.5, 2.5)), 
					truncQ = ifelse(!is.null(options$filtTruncQ), as.numeric(as.character(options$filtTruncQ)), c(0, 0)), 
					truncLen = ifelse(!is.null(options$filtTruncLen), as.numeric(as.character(options$filtTruncLen)), c(0, 0)), 
					trimLeft = ifelse(!is.null(options$filtTrimLeft), as.numeric(as.character(options$filtTrimLeft)), c(0, 0)),
					trimRight = ifelse(!is.null(options$filtTrimRight), as.numeric(as.character(options$filtTrimRight)), c(0, 0)), 
					matchIDs = ifelse(!is.null(options$filtMatchIDs), as.numeric(as.character(options$filtMatchIDs)), FALSE), 
					minLen = ifelse(!is.null(options$filtMinLen), as.numeric(as.character(options$filtMinLen)), c(50, 50)),
					isPaired = TRUE, 
					multithread = ifelse(!is.null(options$multithread), options$multithread, FALSE)) 
		}
		else if ("auto" %in% type & length(type) == 2) {
			temp <- setFastFilter(
						allowOneOff = ifelse(!is.null(options$chimeraDetectionAllowOneOff), options$chimeraDetectionAllowOneOff, FALSE), 
						dadaBandSize = ifelse(!is.null(options$dadaBandSize), as.numeric(as.character(options$dadaBandSize)), 16), 
						dadaOmegaA = ifelse(!is.null(options$dadaOmegaA), as.numeric(as.character(options$dadaOmegaA)), 1e-40), 
						derepN = ifelse(!is.null(options$derepN), as.numeric(as.character(options$derepN)), 1e+06),
						detectionAbundance = ifelse(!is.null(options$chimeraDetectionParentAbundance), as.numeric(as.character(options$chimeraDetectionParentAbundance)), 2), 
						errN = ifelse(!is.null(options$errN), as.numeric(as.character(options$errN)), 1e+08), 
						getChimeraTable = ifelse(!is.null(options$createChimeraDetectionTable), options$createChimeraDetectionTable, FALSE), 
						getErrPDF = ifelse(!is.null(options$saveErrorsPlot), options$saveErrorsPlot, FALSE),
						ignoreNegatives = ifelse(!is.null(options$chimeraDetectionIgnoreNegatives), as.numeric(as.character(options$chimeraDetectionIgnoreNegatives)), 1), 
						inDir = options$pathToData,
						isPaired = TRUE, 
						justConcatenate = ifelse(!is.null(options$mergePairsJustConcatenate), options$mergePairsJustConcatenate, FALSE),
						maxMismatch = ifelse(!is.null(options$mergePairsMaxMismatch), as.numeric(as.character(options$mergePairsMaxMismatch)), 0), 
						maxShift = ifelse(!is.null(options$chimeraDetectionMaxShift), as.numeric(as.character(options$chimeraDetectionMaxShift)), 16),
						mergeSeqs = ifelse(!is.null(options$mergePairs), options$mergePairs, FALSE), 
						minFoldParentOverAbundance = ifelse(!is.null(options$chimeraDetectionMinFoldParentOverabundance), as.numeric(as.character(options$chimeraDetectionMinFoldParentOverabundance)), 1.5), 
						minOverlap = ifelse(!is.null(options$mergePairsMinOverlap), as.numeric(as.character(options$mergePairsMinOverlap)), 12),
						minSampleFraction = ifelse(!is.null(options$chimeraDetectionMinSampleFraction), as.numeric(as.character(options$chimeraDetectionMinSampleFraction)), 0.9), 
						multithread = ifelse(!is.null(options$multithread), options$multithread, FALSE),
						outDir = options$outDir, 
						prefix = ifelse(!is.null(options$projectPrefix), options$projectPrefix, "myproject"),
						returnRejects = ifelse(!is.null(options$mergePairsReturnRejects), options$mergePairsReturnRejects, FALSE), 
						trimOverhang = ifelse(!is.null(options$mergePairsTrimOverhang), options$mergePairsTrimOverhang, FALSE),
						verbose = ifelse(!is.null(options$verbose), options$verbose, FALSE),
						
						maxEE = ifelse(!is.null(options$filtMaxEE), as.numeric(as.character(options$filtMaxEE)), c(2.5, 2.5)), 
						truncQ = ifelse(!is.null(options$filtTruncQ), as.numeric(as.character(options$filtTruncQ)), c(0, 0)), 
						truncLen = ifelse(!is.null(options$filtTruncLen), as.numeric(as.character(options$filtTruncLen)), c(0, 0)), 
						trimLeft = ifelse(!is.null(options$filtTrimLeft), as.numeric(as.character(options$filtTrimLeft)), c(0, 0)),
						trimRight = ifelse(!is.null(options$filtTrimRight), as.numeric(as.character(options$filtTrimRight)), c(0, 0)), 
						matchIDs = ifelse(!is.null(options$filtMatchIDs), as.numeric(as.character(options$filtMatchIDs)), FALSE), 
						minLen = ifelse(!is.null(options$filtMinLen), as.numeric(as.character(options$filtMinLen)), c(50, 50)))
		}
		else {
			stop("Invalid input provided for type parameter")
		}
	}
	else if (type == "auto" & isPaired) {
		# fastPaired object is created
		temp <- setFastPaired(inDir = options$pathToData,
							outDir = options$outDir, 
							mergeSeqs = ifelse(!is.null(options$mergePairs), options$mergePairs, FALSE), 
							trimOverhang = ifelse(!is.null(options$mergePairsTrimOverhang), options$mergePairsTrimOverhang, FALSE),
							minOverlap = ifelse(!is.null(options$mergePairsMinOverlap), as.numeric(as.character(options$mergePairsMinOverlap)), 12),
							maxMismatch = ifelse(!is.null(options$mergePairsMaxMismatch), as.numeric(as.character(options$mergePairsMaxMismatch)), 0), 
							returnRejects = ifelse(!is.null(options$mergePairsReturnRejects), options$mergePairsReturnRejects, FALSE), 
							justConcatenate = ifelse(!is.null(options$mergePairsJustConcatenate), options$mergePairsJustConcatenate, FALSE),
							verbose = ifelse(!is.null(options$verbose), options$verbose, FALSE),
							prefix = ifelse(!is.null(options$projectPrefix), options$projectPrefix, "myproject"),
							isPaired = TRUE, 
							derepN = ifelse(!is.null(options$derepN), as.numeric(as.character(options$derepN)), 1e+06),
							getErrPDF = ifelse(!is.null(options$saveErrorsPlot), options$saveErrorsPlot, FALSE),
							errN = ifelse(!is.null(options$errN), as.numeric(as.character(options$errN)), 1e+08), 
							multithread = ifelse(!is.null(options$multithread), options$multithread, FALSE),
							dadaBandSize = ifelse(!is.null(options$dadaBandSize), as.numeric(as.character(options$dadaBandSize)), 16), 
							dadaOmegaA = ifelse(!is.null(options$dadaOmegaA), as.numeric(as.character(options$dadaOmegaA)), 1e-40), 
							getChimeraTable = ifelse(!is.null(options$createChimeraDetectionTable), options$createChimeraDetectionTable, FALSE), 
							minSampleFraction = ifelse(!is.null(options$chimeraDetectionMinSampleFraction), options$chimeraDetectionMinSampleFraction, 0.9), 
							ignoreNegatives = ifelse(!is.null(options$chimeraDetectionIgnoreNegatives), as.numeric(as.character(options$chimeraDetectionIgnoreNegatives)), 1), 
							minFoldParentOverAbundance = ifelse(!is.null(options$chimeraDetectionMinFoldParentOverabundance), as.numeric(as.character(options$chimeraDetectionMinFoldParentOverabundance)), 1.5), 
							detectionAbundance = ifelse(!is.null(options$chimeraDetectionParentAbundance), as.numeric(as.character(options$chimeraDetectionParentAbundance)), 2), 
							allowOneOff = ifelse(!is.null(options$chimeraDetectionAllowOneOff), options$chimeraDetectionAllowOneOff, FALSE), 
							maxShift = ifelse(!is.null(options$chimeraDetectionMaxShift), as.numeric(as.character(options$chimeraDetectionMaxShift)), 16))
	}
	else if (type == "auto" & !isPaired) {
		temp <- setFastSingle(inDir = options$pathToData,
							outDir = options$outDir, 
							verbose = ifelse(!is.null(options$verbose), options$verbose, FALSE),
							prefix = ifelse(!is.null(options$projectPrefix), options$projectPrefix, "myproject"),
							isPaired = FALSE, 
							derepN = ifelse(!is.null(options$derepN), as.numeric(as.character(options$derepN)), 1e+06),
							getErrPDF = ifelse(!is.null(options$saveErrorsPlot), options$saveErrorsPlot, FALSE),
							errN = ifelse(!is.null(options$errN), as.numeric(as.character(options$errN)), 1e+08), 
							multithread = ifelse(!is.null(options$multithread), options$multithread, FALSE),
							dadaBandSize = ifelse(!is.null(options$dadaBandSize), as.numeric(as.character(options$dadaBandSize)), 16), 
							dadaOmegaA = ifelse(!is.null(options$dadaOmegaA), as.numeric(as.character(options$dadaOmegaA)), 1e-40), 
							getChimeraTable = ifelse(!is.null(options$createChimeraDetectionTable), options$createChimeraDetectionTable, FALSE), 
							minSampleFraction = ifelse(!is.null(options$chimeraDetectionMinSampleFraction), as.numeric(as.character(options$chimeraDetectionMinSampleFraction)), 0.9), 
							ignoreNegatives = ifelse(!is.null(options$chimeraDetectionIgnoreNegatives), as.numeric(as.character(options$chimeraDetectionIgnoreNegatives)), 1), 
							minFoldParentOverAbundance = ifelse(!is.null(options$chimeraDetectionMinFoldParentOverabundance), as.numeric(as.character(options$chimeraDetectionMinFoldParentOverabundance)), 1.5), 
							detectionAbundance = ifelse(!is.null(options$chimeraDetectionParentAbundance), as.numeric(as.character(options$chimeraDetectionParentAbundance)), 2), 
							allowOneOff = ifelse(!is.null(options$chimeraDetectionAllowOneOff), options$chimeraDetectionAllowOneOff, FALSE), 
							maxShift = ifelse(!is.null(options$chimeraDetectionMaxShift), as.numeric(as.character(options$chimeraDetectionMaxShift)), 16))
	}
	else if ((type == "assignTax" | "assignTax" %in% type) & !isPaired) {
		if (length(type) == 1) {
			temp <- setFastAssignTaxa(refDatabase = options$taxDatabase, 
						prefix = ifelse(!is.null(options$projectPrefix), options$projectPrefix, "myproject"),
						minBootstrap = ifelse(!is.null(options$assignTaxMinBootstrap), as.numeric(as.character(options$assignTaxMinBootstrap)), 50), 
						tryComplement = ifelse(!is.null(options$assignTaxTryComplement), options$assignTaxTryComplement, FALSE), 
						showBootstraps = ifelse(!is.null(options$assignTaxOutputBootstraps), options$assignTaxOutputBootstraps, FALSE), 
						taxLevels = ifelse(!is.null(options$assignTaxLevels), options$assignTaxLevels, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")),
						verbose = ifelse(!is.null(options$verbose), options$verbose, FALSE), 
						multithread = ifelse(!is.null(options$multithread), options$multithread, FALSE))
		}
		else if ("auto" %in% type & length(type) == 2) {
			temp <- setFastSingle(inDir = options$pathToData,
								outDir = options$outDir, 
								verbose = ifelse(!is.null(options$verbose), options$verbose, FALSE),
								prefix = ifelse(!is.null(options$projectPrefix), options$projectPrefix, "myproject"),
								isPaired = FALSE, 
								derepN = ifelse(!is.null(options$derepN), as.numeric(as.character(options$derepN)), 1e+06),
								getErrPDF = ifelse(!is.null(options$saveErrorsPlot), options$saveErrorsPlot, FALSE),
								errN = ifelse(!is.null(options$errN), as.numeric(as.character(options$errN)), 1e+08), 
								multithread = ifelse(!is.null(options$multithread), options$multithread, FALSE),
								dadaBandSize = ifelse(!is.null(options$dadaBandSize), as.numeric(as.character(options$dadaBandSize)), 16), 
								dadaOmegaA = ifelse(!is.null(options$dadaOmegaA), as.numeric(as.character(options$dadaOmegaA)), 1e-40), 
								getChimeraTable = ifelse(!is.null(options$createChimeraDetectionTable), options$createChimeraDetectionTable, FALSE), 
								minSampleFraction = ifelse(!is.null(options$chimeraDetectionMinSampleFraction), as.numeric(as.character(options$chimeraDetectionMinSampleFraction)), 0.9), 
								ignoreNegatives = ifelse(!is.null(options$chimeraDetectionIgnoreNegatives), as.numeric(as.character(options$chimeraDetectionIgnoreNegatives)), 1), 
								minFoldParentOverAbundance = ifelse(!is.null(options$chimeraDetectionMinFoldParentOverabundance), as.numeric(as.character(options$chimeraDetectionMinFoldParentOverabundance)), 1.5), 
								detectionAbundance = ifelse(!is.null(options$chimeraDetectionParentAbundance), as.numeric(as.character(options$chimeraDetectionParentAbundance)), 2), 
								allowOneOff = ifelse(!is.null(options$chimeraDetectionAllowOneOff), options$chimeraDetectionAllowOneOff, FALSE), 
								maxShift = ifelse(!is.null(options$chimeraDetectionMaxShift), as.numeric(as.character(options$chimeraDetectionMaxShift)), 16),
								
								refDatabase = options$taxDatabase, 
								minBootstrap = ifelse(!is.null(options$assignTaxMinBootstrap), as.numeric(as.character(options$assignTaxMinBootstrap)), 50), 
								tryComplement = ifelse(!is.null(options$assignTaxTryComplement), options$assignTaxTryComplement, FALSE), 
								showBootstraps = ifelse(!is.null(options$assignTaxOutputBootstraps), options$assignTaxOutputBootstraps, FALSE), 
								taxLevels = ifelse(!is.null(options$assignTaxLevels), options$assignTaxLevels, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")),
								verbose = ifelse(!is.null(options$verbose), options$verbose, FALSE), 
								multithread = ifelse(!is.null(options$multithread), options$multithread, FALSE))
		}
		else {
			stop("Invalid input provided for type parameter")
		}				
	}
	else if ((type == "filter" | "filter" %in% type) & !isPaired) {
		if (length(type) == 1) {
			temp <- setFastFilter(inDir = options$pathToData,
								outDir = options$outDir, 
								filtVerbose = ifelse(!is.null(options$verbose), options$verbose, FALSE),
								prefix = ifelse(!is.null(options$projectPrefix), options$projectPrefix, "myproject"), 
								maxEE = ifelse(!is.null(options$filtMaxEE), as.numeric(as.character(options$filtMaxEE)), c(2.5, 2.5)), 
								truncQ = ifelse(!is.null(options$filtTruncQ), as.numeric(as.character(options$filtTruncQ)), c(0, 0)), 
								truncLen = ifelse(!is.null(options$filtTruncLen), as.numeric(as.character(options$filtTruncLen)), c(0, 0)), 
								trimLeft = ifelse(!is.null(options$filtTrimLeft), as.numeric(as.character(options$filtTrimLeft)), c(0, 0)),
								trimRight = ifelse(!is.null(options$filtTrimRight), as.numeric(as.character(options$filtTrimRight)), c(0, 0)), 
								matchIDs = ifelse(!is.null(options$filtMatchIDs), as.numeric(as.character(options$filtMatchIDs)), FALSE), 
								minLen = ifelse(!is.null(options$filtMinLen), as.numeric(as.character(options$filtMinLen)), c(50, 50)),
								isPaired = TRUE, 
								multithread = ifelse(!is.null(options$multithread), options$multithread, FALSE)) 
		}		
		else if ("auto" %in% type & length(type) == 2) {
			temp <- setFastSingle(inDir = options$pathToData,
								outDir = options$outDir, 
								verbose = ifelse(!is.null(options$verbose), options$verbose, FALSE),
								prefix = ifelse(!is.null(options$projectPrefix), options$projectPrefix, "myproject"),
								isPaired = FALSE, 
								derepN = ifelse(!is.null(options$derepN), as.numeric(as.character(options$derepN)), 1e+06),
								getErrPDF = ifelse(!is.null(options$saveErrorsPlot), options$saveErrorsPlot, FALSE),
								errN = ifelse(!is.null(options$errN), as.numeric(as.character(options$errN)), 1e+08), 
								multithread = ifelse(!is.null(options$multithread), options$multithread, FALSE),
								dadaBandSize = ifelse(!is.null(options$dadaBandSize), as.numeric(as.character(options$dadaBandSize)), 16), 
								dadaOmegaA = ifelse(!is.null(options$dadaOmegaA), as.numeric(as.character(options$dadaOmegaA)), 1e-40), 
								getChimeraTable = ifelse(!is.null(options$createChimeraDetectionTable), options$createChimeraDetectionTable, FALSE), 
								minSampleFraction = ifelse(!is.null(options$chimeraDetectionMinSampleFraction), as.numeric(as.character(options$chimeraDetectionMinSampleFraction)), 0.9), 
								ignoreNegatives = ifelse(!is.null(options$chimeraDetectionIgnoreNegatives), as.numeric(as.character(options$chimeraDetectionIgnoreNegatives)), 1), 
								minFoldParentOverAbundance = ifelse(!is.null(options$chimeraDetectionMinFoldParentOverabundance), as.numeric(as.character(options$chimeraDetectionMinFoldParentOverabundance)), 1.5), 
								detectionAbundance = ifelse(!is.null(options$chimeraDetectionParentAbundance), as.numeric(as.character(options$chimeraDetectionParentAbundance)), 2), 
								allowOneOff = ifelse(!is.null(options$chimeraDetectionAllowOneOff), options$chimeraDetectionAllowOneOff, FALSE), 
								maxShift = ifelse(!is.null(options$chimeraDetectionMaxShift), as.numeric(as.character(options$chimeraDetectionMaxShift)), 16),

								maxEE = ifelse(!is.null(options$filtMaxEE), as.numeric(as.character(options$filtMaxEE)), c(2.5, 2.5)), 
								truncQ = ifelse(!is.null(options$filtTruncQ), as.numeric(as.character(options$filtTruncQ)), c(0, 0)), 
								truncLen = ifelse(!is.null(options$filtTruncLen), as.numeric(as.character(options$filtTruncLen)), c(0, 0)), 
								trimLeft = ifelse(!is.null(options$filtTrimLeft), as.numeric(as.character(options$filtTrimLeft)), c(0, 0)),
								trimRight = ifelse(!is.null(options$filtTrimRight), as.numeric(as.character(options$filtTrimRight)), c(0, 0)), 
								matchIDs = ifelse(!is.null(options$filtMatchIDs), as.numeric(as.character(options$filtMatchIDs)), FALSE), 
								minLen = ifelse(!is.null(options$filtMinLen), as.numeric(as.character(options$filtMinLen)), c(50, 50)))
		}
		else {
			stop("Invalid input provided for type parameter")
		}	
	}
	else if (type == "report") {
		temp <- setFastReport(inDir = options$pathToData,
							outDir = options$outDir, 
							fastqcPath = options$pathToFastqc, 
							installFastqc = ifelse(!is.null(options$installFastqc), options$installFastqc, FALSE),
							numThreads = ifelse(!is.null(options$fastqcThreads), as.numeric(as.character(options$fastqcThreads)), 4), 
							description = ifelse(!is.null(options$fastqcExperimentDescription), options$fastqcExperimentDescription, "My Project"))
	}
	else if (type == "seqdump") {
		temp <- setFastSeqDump(sampleURLs = options$pathToSampleURLs, 
							outDir = options$outDir, 
							sampleList = options$pathToSampleIDs, 
							fastqDumpPath = options$pathToFastqDump)
	
	}
	else if (type == "primertrim") {
		temp <- setFastPrimerTrim(inDir = options$pathToRawFastq,
							outDir = options$pathToNoPrimers,  
							adapterList = options$listOfAdapters)
	}
	else if (type == "qualityplot") {
		temp <- setfastPlotQuality(aggregate = ifelse(!is.null(options$aggregateQual), options$aggregateQual, TRUE), 
							N = ifelse(!is.null(options$qualN), options$qualN, 5e+05))
	}
	else { 
		stop("Error reading config file. Invalid inputs supplied")
	}
	
	# Return object
	return(temp)
}
