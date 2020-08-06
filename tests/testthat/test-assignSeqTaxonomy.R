# R version 3.5.3 (2019-03-11)
# DADA2 version: 1.10.1

library(testthat)
library(dada2) # Install using BiocManager --- installed latticeExtra using devtools::install_version("latticeExtra", version = "0.6-28")

context("Assigning Taxonomy")

source("C:/Users/twuma/FettweisLab/fastq2otu/R/readConfig.R")
source("C:/Users/twuma/FettweisLab/fastq2otu/R/assignSeqTaxonomy.R")

# Provide path to config file
config <- "C:/Users/twuma/FettweisLab/fastq2otu/inst/example-config.yml"
options <- yaml::yaml.load_file(config)

# Update options to test
options$taxDatabase <- "C:/Users/twuma/Downloads/silva_v138_lpsn_accession_matches_v1v4_primersearched_taxa.fa"
options$outDir <- "C:/Users/twuma/Downloads/"

# Custruct sequence table (uses DADA2's built in dataset)
derep1 <- derepFastq(system.file("extdata", "sam1F.fastq.gz", package="dada2"))
derep2 <- derepFastq(system.file("extdata", "sam2F.fastq.gz", package="dada2"))
dada1 <- dada(derep1, tperr1)
dada2 <- dada(derep2, tperr1)
seqtab <- dada2::makeSequenceTable(list(sample1=dada1, sample2=dada2))


# Perform Tests
test_that("Assign taxonomy using custom database", {
  object <- readConfig(config, type = "assignTax")
  
  otutab <- assignSeqTaxonomy(seqtab, object)
  expect_true(!is.null(object))
})
