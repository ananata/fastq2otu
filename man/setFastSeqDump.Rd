% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setup.R
\name{setFastSeqDump}
\alias{setFastSeqDump}
\title{Create a fastSeqDump object
Users have the option of downloading SRA data using fastq-dump or FTP.}
\usage{
setFastSeqDump(
  sampleURLs,
  outDir,
  validate = FALSE,
  sampleList,
  useFastqDump = FALSE,
  fastqDumpPath = NA_character_,
  install_dump = FALSE
)
}
\arguments{
\item{sampleURLs}{OPTIONAL input. Path to text file containing FTP URLs generated from SRA Explorer. If this parameter is provided, then fastqDumpPath is ignored.}

\item{outDir}{Required input. Path to write output files.}

\item{validate}{Default is FALSE. When TRUE, all validation methods are executed}

\item{sampleList}{Required input if fastqDumpPath is specified. New-line delimited text file containin SRA accession IDs.}

\item{useFastqDump}{If TRUE, Fastq-dump will be used in place of wget. Requires fastqDumpPath to be supplied. If not present, fastq-dump present in system.files will be used.}

\item{fastqDumpPath}{Set path to fastq-dump script. If provided, fastq-dump is used in place of wget/curl.

TODO: Add Example}
}
\description{
Create a fastSeqDump object
Users have the option of downloading SRA data using fastq-dump or FTP.
}
