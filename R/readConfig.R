#' Read in a config file (preferred method)
#'
#' This function parses config file to extract data needed to create a FASTQ2OTU class object. It assumes that
#' the file is in YML format and contains default key terms.
#' Execute getSlots("fastq2otu") to get all required arguments or refer to documenation for more details and/or
#' examples.
#'
#' Use when config file is supplied and all required values are provided.
#' @param configFile Path to config file (YML-formatted)
#' @return An S4 object
#' @export
readConfig <- function(configFile) {
	if (!file.exists(configFile)) {
			stop(sprintf("'%s' could not be found, please enter a valid path", configFile))
	}
	options <- yaml::yaml.load_file(configFile)

	# Create new fastq2otu class object based on sequencing method
	if (options$isPaired) {
		temp <- new("fastq2otu_paired", 
					taxDatabase = options$taxDatabase,
					pathToPairedData = options$pathToData,
					projectPrefix = options$projectPrefix,
					outDir = options$outDir,

					pathToSampleIDs = options$pathToSampleIDs,
					runFastqDump = options$runFastqDump,

					trimPrimers = options$trimPrimers,
					listOfAdapters = options$listOfAdapters,
					pathToRawFastq = options$pathToRawFastq,
					pathToNoPrimers = options$pathToNoPrimers,

					filtMaxEE = options$filtMaxEE,
					filtTruncQ = options$filtTruncQ,
					filtTruncLen = options$filtTruncLen,
					filtTrimLeft = options$filtTrimLeft,
					filtTrimRight = options$filtTrimRight,
					filtMultiThread = options$filtMultiThread,
					filtVerbose = options$filtVerbose,
					filtMatchIDs = options$filtMatchIDs,
					filtMinLen = options$filtMinLen,

					doMergeSeqPairs = options$doMergeSeqPairs,
					mergeSeqPairsTrimOverhang = options$mergeSeqPairsTrimOverhang,
					mergeSeqPairsMinOverlap = options$mergeSeqPairsMinOverlap,
					mergeSeqPairsMaxMismatch = options$mergeSeqPairsMaxMismatch,
					mergeSeqPairsReturnRejects = options$mergeSeqPairsReturnRejects,
					mergeSeqPairsJustConcatenate = options$mergeSeqPairsJustConcatenate,
					mergeSeqPairsVerbose = options$mergeSeqPairsVerbose,

					runFastqc = options$runFastqc,
					pathToFastqcResults = options$pathToFastqcResults,
					fastqcThreads = options$fastqcThreads,
					fastqcExperimentDescription = options$fastqcExperimentDescription,
				
					finalSummaryTable = options$finalSummaryTable,

					plotQuality = options$plotQuality,
					qualityPlotPDF = options$qualityPlotPDF,

					derepVerbose = options$derepVerbose,
					derepN = options$derepN,

					errN = options$errN,
					errMultithread = options$errMultithread,
					learnErrorsPlotPDF = options$learnErrorsPlotPDF,

					dadaBandSize = options$dadaBandSize,
					dadaOmegaA = options$dadaOmegaA,

					createChimeraDetectionTable = options$createChimeraDetectionTable,
					chimeraDetectionMinSampleFraction = options$chimeraDetectionMinSampleFraction,
					chimeraDetectionIgnoreNegatives = options$chimeraDetectionIgnoreNegatives,
					chimeraDetectionMinFoldParentOverabundance = options$chimeraDetectionMinFoldParentOverabundance,
					chimeraDetectionParentAbundance = options$chimeraDetectionParentAbundance,
					chimeraDetectionAllowOneOff = options$chimeraDetectionAllowOneOff,
					chimeraDetectionMaxShift = options$chimeraDetectionMaxShift,
					chimeraDetectionMultiThread = options$chimeraDetectionMultiThread,
					chimeraDetectionVerbose = options$chimeraDetectionVerbose,

					assignTaxMinBootstrap = options$assignTaxMinBootstrap,
					assignTaxTryComplement = options$assignTaxTryComplement,
					assignTaxOutputBootstraps = options$assignTaxOutputBootstraps,
					assignTaxLevels = options$assignTaxLevels,
					assignTaxMultiThread = options$assignTaxMultiThread,
					assignTaxVerbose = options$assignTaxVerbose,

					mergeSamples = options$mergeSamples,
					finalMergedTable = options$finalMergedTable)
	} else {
		temp <- new("fastq2otu_single", 
					taxDatabase = options$taxDatabase,
					pathToSingleData = options$pathToSingleData,
					projectPrefix = options$projectPrefix,
					outDir = options$outDir,

					pathToSampleIDs = options$pathToSampleIDs,
					runFastqDump = options$runFastqDump,

					trimPrimers = options$trimPrimers,
					listOfAdapters = options$listOfAdapters,
					pathToRawFastq = options$pathToRawFastq,
					pathToNoPrimers = options$pathToNoPrimers,

					filtMaxEE = options$filtMaxEE,
					filtTruncQ = options$filtTruncQ,
					filtTruncLen = options$filtTruncLen,
					filtTrimLeft = options$filtTrimLeft,
					filtTrimRight = options$filtTrimRight,
					filtMultiThread = options$filtMultiThread,
					filtVerbose = options$filtVerbose,
					filtMatchIDs = options$filtMatchIDs,
					filtMinLen = options$filtMinLen,

					runFastqc = options$runFastqc,
					pathToFastqcResults = options$pathToFastqcResults,
					fastqcThreads = options$fastqcThreads,
					fastqcExperimentDescription = options$fastqcExperimentDescription,
				
					finalSummaryTable = options$finalSummaryTable,

					plotQuality = options$plotQuality,
					qualityPlotPDF = options$qualityPlotPDF,

					derepVerbose = options$derepVerbose,
					derepN = options$derepN,

					errN = options$errN,
					errMultithread = options$errMultithread,
					learnErrorsPlotPDF = options$learnErrorsPlotPDF,

					dadaBandSize = options$dadaBandSize,
					dadaOmegaA = options$dadaOmegaA,

					createChimeraDetectionTable = options$createChimeraDetectionTable,
					chimeraDetectionMinSampleFraction = options$chimeraDetectionMinSampleFraction,
					chimeraDetectionIgnoreNegatives = options$chimeraDetectionIgnoreNegatives,
					chimeraDetectionMinFoldParentOverabundance = options$chimeraDetectionMinFoldParentOverabundance,
					chimeraDetectionParentAbundance = options$chimeraDetectionParentAbundance,
					chimeraDetectionAllowOneOff = options$chimeraDetectionAllowOneOff,
					chimeraDetectionMaxShift = options$chimeraDetectionMaxShift,
					chimeraDetectionMultiThread = options$chimeraDetectionMultiThread,
					chimeraDetectionVerbose = options$chimeraDetectionVerbose,

					assignTaxMinBootstrap = options$assignTaxMinBootstrap,
					assignTaxTryComplement = options$assignTaxTryComplement,
					assignTaxOutputBootstraps = options$assignTaxOutputBootstraps,
					assignTaxLevels = options$assignTaxLevels,
					assignTaxMultiThread = options$assignTaxMultiThread,
					assignTaxVerbose = options$assignTaxVerbose,

					mergeSamples = options$mergeSamples,
					finalMergedTable = options$finalMergedTable)
	}
	# Return object
	return(temp)
}
