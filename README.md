# FASTQ2OTU

FASTQ2OTU was developed as a easy and effective tool for downloading, analyzing, and processing microbiome rRNA gene data obtained from NCBI's SRA database. The package uses many functions from [DADA2](https://github.com/benjjneb/dada2 "Github") to analyze sequence data. The primary objective of FASTQ2OTU is to (1) increase the reproducibility of microbiome analysis and (2) encourage the analysis of archived data to obtain new knowledge. 

FASTQ2OTU's workflow can be broken down into 
#### Advantages of using FASTQ2OTU
* ##### Documentation is simplified
* ##### Integrated workflow
* ##### Outputs are automatically generated
* ##### Easy to use

## Quick Start Guide

#### Directory Overview
```
fastq2otu:.
|   DESCRIPTION
|   fastq2otu.Rproj
|   NAMESPACE
|   README.md
|
+---data
+---exec
|       bbduk.sh
|       fastq-dump
|       retrieve_sra_sequences.sh
|       use_bbduk.sh
|
+---inst
|       bbduk_LISCENSE
|       example-config.yml
|       fastq-dump_LISCENSE
|
\---R
        assignSeqTaxonomy.R
        dadaSeqs.R
        filtTrim.R
        getRowSums.R
        getSeqs.R
        learnSeqErrors.R
        makeSeqsTable.R
        mergeSamples.R
        mergeSeqPairs.R
        plotQuality.R
        readConfig.R
        removeChimeras.R
        removePrimers.R
        runFastqc.R
        saveSeqs.R
        saveTaxonomyTables.R
        setup.R
```
For navigation purposes, the above diagram has been provided as a general schematic of all the files and sub-directories located within the FASTQ2OTU package. 

#### Demo Files
```
# Load package into environoment
library("ananata/fastq2otu")

// Run demo
runDemo()
```
The `runDemo()` function will execute the entire pipeline on a set of a sample data. The example datasets (including config file) can be viewed in the `inst/` directory (TODO), and output files will be written to the user's current working directory. 

## Getting Started
After installing FASTQ2OTU, the following input files and/or directories will be required to begin processing data:
- A YML-formatted config file containing all parameters (more information about the config file can be found below).
- A directory of single or paired-end FASTQ files _OR_ a text file containing SRA ids to download from NCBI. 
- A bit of knowledge about the sequences
    - Ideal trimming parameters
    - Forward and reverse primer lengths
    - If merging, the desired overlap length
    - Working knowledge of [DADA2 workflow](https://benjjneb.github.io/dada2/tutorial.html)

### Prerequisites and Dependencies


### Installation

This application is designed to be lightweight and simple to use. The intended use is via a remote server, however it can also be run using RStudio (the package was written in R 3.5.3) and can be downloaded from Github using [devtools](https://github.com/r-lib/devtools). 

```
# Install devtools 
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("ananata/fastq2otu")

# Load package into environoment
library("ananata/fastq2otu")

```
## Running FASTQ2OTU

### Using a config file

| Variable        | Type           | Default  | Description |
| --------------- |--------------|--------|-------------|
|taxDatabase|Character|N/A| Required. Path to reference taxonomy database. |
|isPaired|Logical| FALSE | Required. Determine whether input sequences are single (FALSE) or paired-end (TRUE) |
|pathToData|Character| N/A | Required. Path to directory containing input FASTQ files. For paired-end data, files containing forward or reverse reads must be in the same directory path. | 
|projectPrefix|Character|myproject| Unique identifier to label output files generated from workflow.|
|outDir|Character|Current working directory|Path to output directory that the contain all output files and documents. |
|runFastqDump|Logical|FALSE|Determines whether fastq-dump should be executed. |
|pathToSampleIDs|Character|N/A|Required if runFastqDump is TRUE. Path to a list of SRA Accession ids to be downloaded using [fastq-dump](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump).|
|trimPrimers|Logical|FALSE|Determines whether bbduk.sh should be executed to trim all adapters present on input FASTQ files |
|listOfAdapters|Character|N/A|Required if trimPrimers is TRUE. Path to a list of adapter sequence to be removed using [bbduk.sh](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/)|
|pathToBBDuk|Character|N/A|Path to bbduk.sh script (should be kept in BBTools root directory. Refer to BBTools [0documentation] (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/) for more information.)|
|pathToRawFastq|Character|N/A|Required if trimPrimers is TRUE. Path to directory containing untrimmed sequences. |
| projectPrefix ||||
| outDir ||||
| pathToData ||||
| isPaired |||| 
| listOfAdapters ||||
| pathToRawFastq ||||
| pathToNoPrimers ||||
| aggregateQual ||||
| qualN||||
|runFastqDump||||
|pathToSampleIDs||||
|runFastqc||||
|pathToFastqcResults||||
|fastqcThreads||||
|fastqcExperimentDescription||||

### Package Walkthrough



## FASTQ2OTU output


## Contributing


## Versioning


## Authors

- **Nana Afia Twumasi-Ankrah** - *Primary developer*


## License

This project is licensed under the GNU GPLv3 License.
This license restricts the usage of this application for non-open sourced systems. Please contact the authors for questions related to relicensing of this software in non-open sourced systems.

## Acknowledgments

We would like to thank the following, without whom this would not have happened:
* Virginia Commonwealth University
* University of Texas - Austin
* TODO: Add all acknowledgements



