library(testthat)
context("Download FASTQ sequences from NCBI SRA database")

root <- "/gpfs_fs/vahmp/users/Fettweis/Atopobium_NTA/Wright_Collaboration/r_package/fastq2otu/"

source(paste0(root, "R/readConfig.R"))
source(paste0(root, "R/setup.R"))
source(paste0(root, "R/getSeqs.R"))


paired_config <- paste0(root, "inst/examples/paired/paired-example_config.yml")
paired_options <- yaml::yaml.load_file(paired_config)

# TODO: Add single-end tests

# TODO: SRA-Explorer website is down (08/30/31) 
#test_that("Use wget to install paired-end data files", {
#  object <- readConfig(paired_config, type = "seqdump")
#  fp <- getSeqs(object)
#  expect_true(fp == paired_options$pathToData)
#})	

test_that("Use fastq-dump to install paired-end data files", {
  object <- readConfig(paired_config, type = c("seqdump")) # TODO: Fix if-else statements so when invalid type combinations are provided a clear(er) error message is given
  fp <- getSeqs(object)
  expect_true(fp == paired_options$pathToData)
})
