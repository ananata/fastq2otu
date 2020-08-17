# R version 3.5.3 (2019-03-11)
# DADA2 version: 1.10.1

library(testthat)
library(dada2) # Install using BiocManager --- installed latticeExtra using devtools::install_version("latticeExtra", version = "0.6-28")

context("Assigning Taxonomy")

# Source files
source("C:/Users/twuma/FettweisLab/fastq2otu/R/readConfig.R")
source("C:/Users/twuma/FettweisLab/fastq2otu/R/assignSeqTaxonomy.R")
source("C:/Users/twuma/FettweisLab/fastq2otu/R/setup.R")

# Provide path to config file (update involved parameters in file)
config <- "C:/Users/twuma/FettweisLab/fastq2otu/inst/example-config.yml"
options <- yaml::yaml.load_file(config)


setFastAssignTaxa()

# Custruct sequence table (uses DADA2's built in dataset)
derep1 <- derepFastq(system.file("extdata", "sam1F.fastq.gz", package="dada2"))
derep2 <- derepFastq(system.file("extdata", "sam2F.fastq.gz", package="dada2"))

dada1 <- dada(derep1, tperr1)
dada2 <- dada(derep2, tperr1)
test.seqtab <- dada2::makeSequenceTable(list(sample1=dada1, sample2=dada2))


# Perform Tests
test_that("Create fastAssignTaxa object using config file", {
  object <- readConfig(config, type = "assignTax")
  expect_true(!is.null(object))
})

test_that("Assign taxonomy using custom database", {
  otuTab <- assignSeqTaxonomy(seqtab = test.seqtab, object = readConfig(config, type = "assignTax"))
  expect_true(length(colnames(otuTab) == length(options$assignTaxLevels)))
})


