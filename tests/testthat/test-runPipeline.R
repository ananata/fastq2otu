library(testthat)
context("Run entire pipeline")

root <- "/fastq2otu/"

source(paste0(root, "R/readConfig.R"))
source(paste0(root, "R/getSeqs.R"))
source(paste0(root, "R/assignSeqTaxonomy.R"))
source(paste0(root, "R/dadaSeqs.R"))
source(paste0(root, "R/filtTrim.R"))
source(paste0(root, "R/getRowSums.R"))
source(paste0(root, "R/learnSeqErrors.R"))
source(paste0(root, "R/makeSeqsTable.R"))
source(paste0(root, "R/mergeSamples.R"))
source(paste0(root, "R/mergeSeqPairs.R"))
source(paste0(root, "R/plotQuality.R"))
source(paste0(root, "R/removeChimeras.R"))
source(paste0(root, "R/runFastqc.R"))
source(paste0(root, "R/saveSeqs.R"))
source(paste0(root, "R/saveTaxonomyTables.R"))
source(paste0(root, "R/setup.R"))
source(paste0(root, "R/trimAdapters.R"))
source(paste0(root, "R/runPipeline.R"))


paired_config <- paste0(root, "inst/examples/paired/my_paired-example_config.yml")
paired_options <- yaml::yaml.load_file(paired_config)

message("Check: ", length(as.numeric(as.character(paired_options$filtMaxEE))))

# First Run: Stopped at syntax error in plotQuality.R script
# getSeqs.R executed in full, but prompt to create output directory did not work (had to manually create output directory in order for function to execute)
# runPipeline(configFile = paired_config, isPaired = TRUE, plotQuality = TRUE, mergeSamples = TRUE, downloadSeqs = TRUE, trimAdapters = FALSE, generateReport = FALSE)

# Second Run: Sequences were already downloaded in first run (downloadSeqs is now set to FALSE)
# runPipeline parameter names were changed to differentiate from function names
# getSeqs.R executed in full despite getDownloadedSeqs being FALSE
# Error message: Error in packageVersion("dada2") : there is no package called ‘dada2’
runPipeline(configFile = paired_config, isPaired = TRUE, getQuality = TRUE, getMergedSamples = TRUE, getDownloadedSeqs = FALSE,
                                getTrimmedAdapters = FALSE, getGeneratedReport = FALSE)


# Create object - simplifies debugging process
object <- readConfig(paired_config, isPaired = TRUE, type = c('auto', 'filter'))
object
