# FASTQ2OTU

FASTQ2OTU was developed as a easy and effective tool for downloading, analyzing, and processing microbiome rRNA gene data obtained from NCBI's SRA database. The package uses many functions from [DADA2](https://github.com/benjjneb/dada2 "Github") to analyze sequence data. The primary objective of FASTQ2OTU is to (1) increase the reproducibility of microbiome analysis and (2) encourage the analysis of archived data to obtain new knowledge. 


#### Advantages of using FASTQ2OTU
* ###### Documentation is simplified
* ###### Outputs are automatically generated
* ###### Easy to use

## Quick Start Guide

#### Example
```
# Load package into environoment
library("ananata/fastq2otu")

// TODO: Add sample data as extdata
configFile <- "path/to/config.yml"

# Create custom object with file
myParams <- readConfig(configFile)

# Execute project workflow
run_fastq2otu(myParams)

```

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

This application is designed to be lightweight and simple to use.  The intended use is via a remote server, however it can also be run using RStudio (the package was written in R 3.5.3) and can be downloaded from Github using [devtools](https://github.com/r-lib/devtools). 

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

The user can set several parameters using environment variables passed into the container at runtime. The environment variables that can be passed are as follows:

| Variable        | Type           | Default  | Description |
| --------------- |:--------------:|:--------:|-------------|

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

