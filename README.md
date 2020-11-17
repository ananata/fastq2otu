# FASTQ2OTU
> :warning: **This package is still under active development**

FASTQ2OTU was developed as a easy and effective tool for downloading, analyzing, and processing large-scale microbiome rRNA gene data obtained from NCBI's SRA database. The package uses many functions from [DADA2](https://github.com/benjjneb/dada2 "Github") to analyze sequence data. The primary objective of FASTQ2OTU is to (1) increase the reproducibility of microbiome analysis and (2) encourage the analysis of archived data to obtain new knowledge. This 

FASTQ2OTU's workflow can be broken down into multiple stages:
1. Get Sequences - Sequences can be downloaded using [FASTQDUMP](https://ncbi.github.io/sra-tools/fastq-dump.html) or `wget` (if FTP links are available).
2. Plot Quality Distribution - Generate a figure that shows the quality distribution of the dataset. This step can be run independently if necessary. 
3. Filter and Trim
4. Learn Errors and Denoise
5. Find and Remove Chimeric Sequences
6. Merge Paired-End Sequences
7. Assign Taxonomy
8. Merge OTU Tables - Merges individual OTU tables to create a single table that can used in downstream analyses. 

## Advantages of using FASTQ2OTU
* Documentation is automated (all inputs and outputs are recorded)
* Integrated workflow
* Outputs are automatically generated
* Easy to use

## Directory Overview
```
fastq2otu:.
|   DESCRIPTION
|   fastq2otu.Rproj
|   NAMESPACE
|   README.md
|
+---inst
|   +---bash
|   |       reformat_fastq.sh
|   |       retrieve_sra_sequences.sh
|   |       use_bbduk.sh
|   |
|   \---examples
|       |   
|       |
|       +---paired
|       |       paired-example_config.yml
|       |       paired-example_SRR_Acc_List.txt
|       |       paired-example_SRR_Url_List.txt     
|       |
|       \---single
|               single-example_config.yml
|               single-example_SRR_ACC_List.txt
|
+---man
|       assignSeqTaxonomy.Rd
|       check_assign_tax.Rd
|       check_fastPaired.Rd
|       check_fastq2otu.Rd
|       check_fastSingle.Rd
|       check_filt_params.Rd
|       check_primer_trim.Rd
|       check_seq_dump.Rd
|       dadaSeqs.Rd
|       fastAssignTaxa-class.Rd
|       fastFilter-class.Rd
|       fastPaired-class.Rd
|       fastPlotQuality-class.Rd
|       fastPrimerTrim-class.Rd
|       fastq2otu-class.Rd
|       fastReport-class.Rd
|       fastSeqDump-class.Rd
|       fastSingle-class.Rd
|       filtTrim.Rd
|       getRowSums.Rd
|       getSeqs.Rd
|       makeSeqsTable.Rd
|       mergeSamples.Rd
|       plot_quality.Rd
|       readConfig.Rd
|       removeChimeras.Rd
|       runFastqc.Rd
|       runPipeline.Rd
|       saveSeqs.Rd
|       saveTaxonomyTables.Rd
|       setFastAssignTaxa.Rd
|       setFastFilter.Rd
|       setFastPaired.Rd
|       setfastPlotQuality.Rd
|       setFastPrimerTrim.Rd
|       setFastReport.Rd
|       setFastSeqDump.Rd
|       setFastSingle.Rd
|       trimAdapters.Rd
|
+---R
|       assignSeqTaxonomy.R
|       filtTrim.R
|       getRowSums.R
|       getSeqs.R
|       mergeSamples.R
|       mergeSeqPairs.R
|       plotQuality.R
|       readConfig.R
|       removeChimeras.R
|       runFastqc.R
|       runPipeline.R
|       saveSeqs.R
|       saveTaxonomyTables.R
|       setup.R
|       trimAdapters.R
|
\---tests
    \---testthat
            test-assignSeqTaxonomy.R
            test-dadaSeqs.R
            test-filtTrim.R
            test-getRowSums.R
            test-getSeqs.R
            test-mergeSamples.R
            test-mergeSeqPairs.R
            test-plotQuality.R
            test-readConfig.R
            test-removeChimeras.R
            test-runFastqc.R
            test-runPipeline.R
            test-saveSeqs.R
            test-saveTaxonomyTables.R
            testthat.R
```
For navigation purposes, the above diagram has been provided as a general schematic of all the files and sub-directories located within the FASTQ2OTU package. 

## Getting Started
After installing FASTQ2OTU, the following input files and/or directories will be required to begin processing data:
- A YML-formatted config file containing all parameters (more information about the config file can be found below).
- A directory of single or paired-end FASTQ files _OR_ a text file containing SRA ids to download from NCBI. 
- A bit of knowledge about the sequences
    - Ideal trimming parameters
    - Forward and reverse primer lengths
    - If merging, the desired overlap length
    - Working knowledge of [DADA2 workflow](https://benjjneb.github.io/dada2/tutorial.html)

## Execute pipeline
```
# Load package into environoment
library("ananata/fastq2otu")

# Path to config file
paired_config <- "path/to/my_paired-example_config.yml"

# Run pipeline
runPipeline(configFile = paired_config, isPaired = TRUE, getQuality = TRUE, getMergedSamples = TRUE, getDownloadedSeqs = TRUE, getGeneratedReport = FALSE)
```
The `runPipeline()` function will allow the the entire DADA2 pipeline to be run. The parameters in the function allows users to specify which steps of the pipeline they would like to execute. The following table provides a description of each parameter and the action(s) it controls.  

| Parameter       | Description | Directions |
| --------------- |-------------|------------|
|configFile | Path to YML-file containing all user inputs | The file must be formatted with the correct variable names (please refer to template)|
|isPaired | TRUE if handling paired-end data and FALSE if handling single-end. | Please note that paired-end and single-end data must be processed seperately (the package cannot analyze both datatypes simultaneously). |
| getQuality | TRUE if you would like to generate a quality distribution plot and FALSE if you would like to skip the step. | This step can be run independently. |
| getMergedSamples | TRUE if you would like to generate a merged sample table and FALSE if you would like to skip the step. | Generates a single table containing data from all samples. |
| getDownloadedSeqs | TRUE if you would like to use `fastq-dump` or `wget` to download data directly from NCBI's SRA database. | Requires a text file containing all SRA sample IDs or FTP download links |
| getGeneratedReport | If TRUE, a FASTQC report is generate using the FASTQCR R-package | This step can also be run independently. |

### Quick Start Guide
DADA2 is an R package that allows users to preform high-resolution taxanomy analyses from FASTQ files. This package will allow most users to analyse datasets using the DADA2 pipeline. This procedure will cover some basics of R programming, installing and running the package on R server, and interpreting some of the outputs generated.
There are two objectives for this document:
1. Introduce new users to DADA2’s functions;
2. To set-up a pipeline for 16S rRNA analyses of target bacterial isolates.

### Plot Quality Distribution
DADA2’s `plotQualityProfile()` function creates a plot(s) that visualizes the overall distribution of quality scores within a dataset. Users can use the plots to make informed decisions about how they would like their data to be processed (i.e. filtering and trimming). The generated quality graphs show colored lines that signify different statistics. 
* Green is the mean quality score for all reads in a single dataset 
* Orange is the median 
* Dashed orange lines demarcate the 25th and 75th quantiles.

### Merging Samples
Sequence tables generated by DADA2’s `makeSequenceTable()` function are formatted as single-row matrices (contain only one row), with consensus sequences as column headings and read counts as elements in the row. OTU Tables (given by DADA2's `assignTaxonomy()` or `assignSpecies()` function) contain taxonomic assignments and sequence variants (ASV). FASTQ2OTU's `mergeSamples` function will merge data from sequence and OTU tables obtained from different samples to generate a single table. The final table can be used to make inter-sample comparisons that may inform downstream analyses.

#### Downloading Data from NCBI
Public data can be accessed from NCBI’s [SRA website](https://www.ncbi.nlm.nih.gov/sra) . To view datasets, enter a project ID (i.e. PRJEB8073), click "Search" and select “Send results to Run Selector" link to view the results interactively. To access the Run Selector tool directly, the following [link](https://trace.ncbi.nlm.nih.gov/Traces/study/?go=home) can also be used. To obtain a list of all SRA accession IDs within a given project, click the "Accession List" button in middle the "Select" panel and wait for the text file to be downloaded.

##### Using FASTQ-DUMP
To download datasets using NCBI's fastq-dump utility, download the [sra-toolkit](https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/) from NCBI and obtain the path to the `fastq-dump` tool. Record the paths to the `fastq-dump` script the SRA accession list in the config file (described below). Make sure to set the `getDownloadedSeqs` parameter to TRUE, when executing the `runPipeline()` function. 

##### Using WGET
To download datasets from SRA using `wget`, navigate to [SRA-Explorer](https://sra-explorer.info/) and input your project ID. Once you click the search icon, a table should appear at the bottom of the window. Select all rows in the tables and store the results by clicking the blue "Add to Collection" button on the right. Please not that the search only outputs a certain number of results each time (with the max being 500). In order to obtain data on more than 500 samples, you must update the "Start at Record" text box after each search. Once you have stored all your samples in your collection, click the shopping cart icon on the top right. Click the tab that says "Raw FastQ Download URLs" and select the download link. Record the path to the downloaded text file in the config file and make sure to set the `getDownloadedSeqs` parameter to TRUE, when executing the `runPipeline()` function. 

#### Generating FASTQCR Report
When `getGeneratedReport` is TRUE, [FASTQCR](https://cran.r-project.org/web/packages/fastqcr/readme/README.html) is used to generate a FASTQCR report that provided quality  information about all input data. 

## Installation
This application is designed to be lightweight and simple to use. The intended use is via a remote server, however it can also be run using RStudio (the package was written in R 3.5.3) and can be downloaded from Github using [devtools](https://github.com/r-lib/devtools). 

```
# Install devtools 
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("ananata/fastq2otu")

# Load package into environoment
library("ananata/fastq2otu")

```

### Using a config file

| Variable        | Type         | Default  | Description |
| --------------- |--------------|----------|-------------|
|projectPrefix |Character|"myproject"| Prefix to append to newly created files (i.e. <myproject>_filtered_files/ is created to store filtered files)|
|outDir|Character|Current working directory|Path to output directory that the contain all output files and documents.|
|pathToData |Character|N/A|Path to directory storing all input data.|
|verbose | Logical| FALSE | Sets `verbose` parameter for all functions |
|multithread | Logical| FALSE | Sets the `multithread` parameter for all functions |
|pathToSampleIDs|Character|N/A|The path to a text file containing SRA Accession IDs.|
|fastaPattern | Character| ^.*[1,2]?.fastq(.gz)?$ | Regex pattern to use when parsing directories for FASTQ files. |
|aggregateQual |Logical|N/A|Provide TRUE if you would like to aggregate your quality profile diagram. |
|qualN||Numeric|0|Enter the number of bases to sample to learn seqence error rates.|
|useFastqDump|Logical|FALSE|Provide TRUE if you would like to download sequences using a locally installed version of SRA's FASTQDUMP|
|pathToFastqDump|Character|N/A| Path to fastq-dump script. Required if useFastqDump parameter is TRUE. |
|pathToSampleURLs|Character|N/A| Path to text file containing FTP download links. |
|pathToFastqc|Character|N/A| Path to fastqc software. Required to use FASTQCR|
|installFastqc|Logical|FALSE|If TRUE, FASTQC will be automatically downloaded into the users home directory. Unless an input for pathToFastqc is provided, then the new download will overwrite the older version. |
|pathToFastqcResults|Character|Path to the directory storing the FASTQC reports. |
|taxDatabase|Character|N/A| Required. Path to reference taxonomy database. |

## Authors

- **Nana Afia Twumasi-Ankrah** - *Primary developer*

## License

This project is licensed under the GNU GPLv3 License.
This license restricts the usage of this application for non-open sourced systems. Please contact the authors for questions related to relicensing of this software in non-open sourced systems.

## Acknowledgments

We would like to thank the following:
* DADA2's Original Developers (Callahan Lab)
* Virginia Commonwealth University - Vaginal Microbiome Consortium
* University of Texas - Austin
* Emory University




