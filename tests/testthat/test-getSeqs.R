library(testthat)
context("Download FASTQ sequences from NCBI SRA database")

source("readConfig.R")
source("setup.R")
source("getSeqs.R")


config <- #TODO
options <- yaml::yaml.load_file(config)


test_that("Use wget to install files", {
  object <- readConfig(config, type = "seqdump")
  fp <- getSeqs(object)
  expect_true(fp == options$pathToData)
})

test_that("Use fastq-dump to install files", {
  object <- readConfig(config, type = c("auto", "seqdump"))
  fp <- getSeqs(object)
  expect_true(fp == options$pathToData)
})
