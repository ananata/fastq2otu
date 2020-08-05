library(testthat)
context("Parse Config File")

source("C:/Users/twuma/FettweisLab/fastq2otu/R/readConfig.R")
source("C:/Users/twuma/FettweisLab/fastq2otu/R/setup.R")


config <- "C:/Users/twuma/FettweisLab/fastq2otu/inst/example-config.yml"
options <- yaml::yaml.load_file(config)

test_that("Create fastSingle object", {
  object <- readConfig(config, isPaired = FALSE, type = "auto")
  expect_true(!is.null(object))
})

test_that("Create fastDouble object", {
  object <- readConfig(config, isPaired = TRUE, type = "auto")
  expect_true(!is.null(object))
})

test_that("Create fastAssignTaxa object", {
  object <- readConfig(config, type = "assignTax")
  expect_true(!is.null(object))
})

test_that("Create fastReport object", {
  object <- readConfig(config, type = "report")
  expect_true(!is.null(object))
})

test_that("Create fastFilter object", {
  object <- readConfig(config, type = "filter")
  expect_true(!is.null(object))
})

test_that("Create fastSeqDump object", {
  object <- readConfig(config, type = "seqdump")
  expect_true(!is.null(object))
})

test_that("Create fastPrimerTrim object", {
  object <- readConfig(config, type = "primertrim")
  expect_true(!is.null(object))
})

test_that("Create fastPlotQuality object", {
  object <- readConfig(config, type = "qualityplot")
  expect_true(!is.null(object))
})




