library(testthat)
context("Parse Config File")

source("readConfig.R")
source("setup.R")


config <- #TODO
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

test_that("Create fastSingle object with Filtering parameters", {
  object <- readConfig(config, type = c('auto', 'filter'))
  expect_true(!is.null(object))
})

test_that("Create fastSingle object with Assign Taxonomy parameters", {
  object <- readConfig(config, type = c('auto', 'assignTax'))
  expect_true(!is.null(object))
})

# Flipping order of list
test_that("Create fastSingle object with Filtering parameters", {
  object <- readConfig(config, type = c('filter', 'auto'))
  expect_true(!is.null(object))
})

test_that("Create fastSingle object with Assign Taxonomy parameters", {
  object <- readConfig(config, type = c('assignTax', 'auto'))
  expect_true(!is.null(object))
})




