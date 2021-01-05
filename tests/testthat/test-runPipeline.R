library(testthat)
context("Run entire pipeline")

root <- "/gpfs_fs/vahmp/users/Fettweis/Atopobium_NTA/Wright_Collaboration/r_package/fastq2otu/"

source(paste0(root, "R/readConfig.R"))
source(paste0(root, "R/getSeqs.R"))
source(paste0(root, "R/assignSeqTaxonomy.R"))
source(paste0(root, "R/filtTrim.R"))
source(paste0(root, "R/getRowSums.R"))
#source(paste0(root, "R/makeSeqsTable.R"))
source(paste0(root, "R/mergeSamples.R"))
source(paste0(root, "R/mergeSeqPairs.R"))
source(paste0(root, "R/plotQuality.R"))
source(paste0(root, "R/removeChimeras.R"))
source(paste0(root, "R/runFastqc.R"))
source(paste0(root, "R/saveSeqs.R"))
source(paste0(root, "R/saveTaxonomyTables.R"))
source(paste0(root, "R/setup.R"))
source(paste0(root, "R/trimAdapters.R"))
source(paste0(root, "R/modified_runPipeline.R"))


paired_config <- paste0(root, "inst/examples/paired/my_paired-example_config.yml")
paired_options <- yaml::yaml.load_file(paired_config)

message("Running Paired Analysis")
runPipeline(configFile = paired_config, isPaired = TRUE, getQuality = TRUE, getMergedSamples = TRUE, getDownloadedSeqs = FALSE, 
	getTrimmedAdapters = FALSE)

single_config <- paste0(root, "inst/examples/single/my_single-example_config.yml")
single_options <- yaml::yaml.load_file(single_config)

message("Running Single-End Analysis")
runPipeline(configFile = single_config, isPaired = FALSE, getQuality = TRUE, getMergedSamples = TRUE, getDownloadedSeqs = FALSE,
	getTrimmedAdapters = FALSE)

