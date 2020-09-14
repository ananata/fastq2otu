library(testthat)
context("Download FASTQ sequences from NCBI SRA database")

root <- "/fastq2otu/"

source(paste0(root, "R/readConfig.R"))
source(paste0(root, "R/setup.R"))
source(paste0(root, "R/getSeqs.R"))


paired_config <- paste0(root, "inst/examples/paired/my_paired-example_config.yml")
paired_options <- yaml::yaml.load_file(paired_config)

# TODO: Add single-end tests

# TODO: SRA-Explorer website is down (08/30/31) 
# Update: Posted example output from SRA-Explorer (see subdirectories of fastq2otu/inst/examples/)

test_that("Use wget to install paired-end data files", {
  object <- readConfig(paired_config, type = "seqdump")
  fp <- getSeqs(object, useFastqDump = FALSE)
  expect_true(fp == paired_options$pathToData)
})	

test_that("Use fastq-dump to install paired-end data files", {
  object <- readConfig(paired_config, type = c("seqdump"))
  fp <- getSeqs(object, useFastqDump = TRUE)
  expect_true(fp == paired_options$pathToData)
})
