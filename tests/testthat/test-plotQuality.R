library(testthat)
context("Generate quality distribution plots")

root <- "/fastq2otu/"

source(paste0(root, "R/readConfig.R"))
source(paste0(root, "R/setup.R"))
source(paste0(root, "R/plotQuality.R"))


paired_config <- paste0(root, "inst/examples/paired/my_paired-example_config.yml")
paired_options <- yaml::yaml.load_file(paired_config)

# TODO: Add single-end tests

test_that("Generate aggregated quality plots", {
  object <- readConfig(paired_config, type = "seqdump")
  fp <- getSeqs(object, useFastqDump = FALSE)
  expect_true(fp == paired_options$pathToData)
})

test_that("Generate non-aggregated quality plots", {
  object <- readConfig(paired_config, type = c("seqdump"))
  fp <- getSeqs(object, useFastqDump = TRUE)
  expect_true(fp == paired_options$pathToData)
})

