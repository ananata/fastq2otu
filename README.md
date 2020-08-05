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
        run_fastq2otu.R
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
devtools::install("ananata/fastq2otu")

# Load package into environoment
library("ananata/fastq2otu")

```
## Running FASTQ2OTU

#### Using a config file

| Variable        | Type           | Default  | Description |
| --------------- |--------------:|:--------:|-------------|
| projectPrefix ||||
| outDir ||||
| pathToData ||||
| isPaired |||| 
| listOfAdapters ||||
| pathToRawFastq ||||
| pathToNoPrimers ||||
| aggregateQual ||||
| qualN||||
|runFastqDump: false
pathToSampleIDs: path/to/SRR_Acc_List.txt
runFastqc: true
pathToFastqcResults: path/to/fastqc_results
fastqcThreads: 4
fastqcExperimentDescription: "16S rRNA amplicon data"

# === Create summary table that displays changes to total read counts ===
finalSummaryTable: filt_summary_table.csv

# === Dereplicate reads to keep only unique sequences ===
derepVerbose: false
derepN: 1e+06

# === Detect and learn error patterns in sequences ===
errN: 1e+08
errMultithread: false
saveErrorsPlot: false

# === Denoise data ===
dadaBandSize: 16
dadaOmegaA: 1e-40

# === Find chimeric sequences ===
createChimeraDetectionTable: false
chimeraDetectionMinSampleFraction: 0.9
chimeraDetectionIgnoreNegatives: 1
chimeraDetectionMinFoldParentOverabundance: 1.5
chimeraDetectionParentAbundance: 2
chimeraDetectionAllowOneOff: false
chimeraDetectionMaxShift: 16
chimeraDetectionMultiThread: false
chimeraDetectionVerbose: false

# === Filtering (Single-end data example)
filtMaxEE: 2.5
filtTruncQ: 0
filtTruncLen: 0
filtTrimLeft: 0
filtTrimRight: 0
filtMultiThread: true
filtVerbose: true
filtMatchIDs: false
filtMinLen: 50

# === Filtering (Paired-end data example)
filtMaxEE: [2.5, 2.5]
filtTruncQ: [0, 0]
filtTruncLen: [0, 0]
filtTrimLeft: [0, 0]
filtTrimRight: [0, 0]
filtMultiThread: true
filtVerbose: true
filtMatchIDs: false
filtMinLen: [50, 50]

# === Merge Paired-end Reads ===
mergePairs: false
mergePairsTrimOverhang: true
mergePairsMinOverlap: 12
mergePairsMaxMismatch: 0
mergePairsReturnRejects: false
mergePairsJustConcatenate: false
mergePairsVerbose: false

# === Assign Taxonomy ===
taxDatabase: path/to/ref
assignTaxMinBootstrap: 50
assignTaxTryComplement: false
assignTaxOutputBootstraps: false
assignTaxLevels: ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
assignTaxMultiThread: true
assignTaxVerbose: true

#### Package Walkthrough



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

---------------------------------------------------------------------------------------------------------------------

