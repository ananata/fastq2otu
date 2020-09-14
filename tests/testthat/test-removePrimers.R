context("Use bbduk.sh to trim adapters from sequences")

root <- "/fastq2otu/"

source(paste0(root, "R/readConfig.R"))
source(paste0(root, "R/setup.R"))
source(paste0(root, "R/removePrimers.R"))


paired_config <- paste0(root, "inst/examples/paired/my_paired-example_config.yml")
paired_options <- yaml::yaml.load_file(paired_config)

# TODO: Add single-end tests
# TOD0: Test using bbduk.sh path

test_that("Remove primers from paired-end data", {
  object <- readConfig(paired_config, type = "primertrim")
  fp <- removePrimers(object)
  expect_true(fp == paired_options$pathToNoPrimers)
})

test_that("Remove primers from single-end data", {
  object <- readConfig(paired_config, type = c("primertrim"))
  fp <- removePrimers(object)
  expect_true(fp == paired_options$pathToNoPrimers)
})

