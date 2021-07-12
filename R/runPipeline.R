#' Import all methods
#' @import methods
NULL


#' runPipeline
#'
#' The analysis takes place in multiple steps beginning at the creation of a central output directory. Once the directory
#' is created, a log-file is initilized that will contain all messages produced by DADA2. Accompanying the log file is a summary
#' table documenting changes in read frequency following filtering and trimming procedures. This summary table can be used
#' as a reference when determining the best way to modify the parameters entered into the different functions.
#'
#' @param plotQuality Default is TRUE. Sequence quality distribution is plotted using DADA2's plotQualityProfile method.
#' @param isPaired Default if FALSE. If TRUE workflow for paired end data is executed.
#' @param getMergeSamples Default is TRUE. If FALSE, sample OTU tables will not be merged across.
#' @param configFile Path to config file (YML-formatted).
#' @param downloadSeqs Default is FALSE. If TRUE, users will be able to retrieve FASTQ files from SRA using the fastq2otu's getSeqs method.
#' @param trimAdapters Default is FALSE. If TRUE, users will be able to remove adapters sequences from data using BBTools bbduk.sh script.
#' @importFrom yaml yaml.load_file
#' @import dada2
#'
#' @export

# ======
# MAIN METHOD
# ======
runPipeline <- function(configFile, isPaired = FALSE, getQuality = TRUE, getMergedSamples = TRUE, getDownloadedSeqs = FALSE,
		getTrimmedAdapters = FALSE) {
  if (!file.exists(configFile)) {
    stop(sprintf("'%s' could not be found, please enter a valid path", configFile))
  }
  
  # Parse YAML file
  options <- yaml::yaml.load_file(configFile)
  
  # Get path to input directory (if valid)
  fp <- options$pathToData
  
  # Get and set path to output directory
  out <- options$outDir
  
  # Get current working directory
  curr.wd <- getwd()
  
  # Change working directory to output directory
  path <- out
  if (dir.exists(path)) {
    base::setwd(path)
  } else {
    #stop("Could not find the following path: ", path)

    # Create directory if it does not exist
    if (!dir.exists(path)) {
        response <- readline(prompt=paste0("Are you sure you want to create ", path, "? <y/N> : "))
        if (response %in% c("Yes", "Y", "y", "yes")) {
          dir.create(path)
          message("Created ", path)
        } else {
          stop("Program was stopped. Could not create directory")
        }
    }

  }

  # Save all outputs to file
  if (!is.null(options$projectPrefix)) {
    log.file <- file.path(out, paste0(options$projectPrefix, "_fastq2otu_output.log"))
  } else {
    log.file <- file.path(out, "fastq2otu_output.log")
    options$projectPrefix <- "fastq2otu_project" # Changed to default
  }
    
  # Print the date
  write(paste0("Date: ", Sys.Date()), log.file)
  write(paste0("Analyzing data for ", options$projectPrefix), log.file, append = TRUE)
  
  # Print package versions
  write(paste0("R version: ", version$version.string), log.file, append = TRUE)
  write(paste0("DADA2 Version: ", packageVersion("dada2")), log.file, append = TRUE)
  write(paste0("YAML Version: ", packageVersion("yaml")), log.file, append = TRUE)
  
  if (getDownloadedSeqs == TRUE) {
    write("==== Downloading sequences from NCBI ====", log.file, append = TRUE)
    message("==== Downloading sequences from NCBI ====")
    object <- readConfig(configFile, type = "seqdump")
    fp <- getSeqs(object)
  }

  # Remove primers and update path
  if (getTrimmedAdapters) {
    write("==== Removing Primers ====", log.file, append = TRUE)
    message("==== Removing Primers ====")
    object <- readConfig(configFile, type = "primertrim")
    fp <- trimAdapters(object)
  } else {
    fp <- options$pathToData
  }
    
  # Plot Quality Distribution and save object
  if (getQuality) {
    write("==== Plotting quality distribution BEFORE trimming ====", log.file, append = TRUE)
    message("==== Plotting quality distribution BEFORE trimming ====")
    plot.object <- readConfig(configFile, type = "qualityplot")
    if (is.null(plot.object)) {
      message("Unable to generate quality plot object")
      stop("Unable to generate quality plot object")
    }

    plots <- plotQuality(fp, options$projectPrefix, plot.object)
    
    # Write PDF files in output directory
    if (!is.null(plots)) {
      write(paste0("Created: ", paste0(options$projectPrefix, "_fastq2otu_quality_plots_BEFORE.pdf")), log.file, append = TRUE)
      message("Created: ", paste0(options$projectPrefix, "_fastq2otu_quality_plots_BEFORE.pdf"))
      pdf(file = file.path(out, paste0(options$projectPrefix, "_fastq2otu_quality_plots_BEFORE.pdf")))
	      print(plots)
      dev.off()
    } else {
      stop("Unable to write quality plots to PDF") # Should rarely execute. Mostly for debugging purposes
    }
  }
  
  if (isPaired) {
    # Get file extension pattern
    if (!is.null(options$fastaPattern) & length(options$fastaPattern) == 2) {
      REGEX_PAT <- options$fastaPattern
    } else {
      message("Missing Input Warning: fastaPattern defaults used.")
      REGEX_PAT <- c("*_1.fastq(.gz)?", "*_2.fastq(.gz)?")
    }

    # === Analyze as paired-end data ===
    sample.names <- na.omit(sapply(strsplit(basename(sort(list.files(fp, pattern = REGEX_PAT[1], full.names = TRUE))),REGEX_PAT[1]), `[`, 1))
    amplicons <- paired_analysis(fp = fp, sample.names = sample.names, file = configFile, getQuality = getQuality, REGEX_1 = REGEX_PAT[1], REGEX_2 = REGEX_PAT[2], logFile = log.file)

  } else {
    # Get file extension pattern
    if (!is.null(options$fastaPattern) & length(options$fastaPattern) == 1) {
      REGEX_PAT <- options$fastaPattern
    } else {
      REGEX_PAT <- "*.fastq(.gz)?$"
    }

    # === Analyze single-end data ===
    sample.names <- na.omit(sapply(strsplit(basename(sort(list.files(fp, pattern = REGEX_PAT, full.names = TRUE))), REGEX_PAT), `[`, 1))
    amplicons <- single_analysis(fp = fp, sample.names = sample.names, file = configFile, getQuality = getQuality, REGEX = REGEX_PAT, logFile = log.file)

  }
  
  # Merge tables
  if (getMergedSamples) {
    write("==== Merging OTU and Sequence Tables ====", log.file, append = TRUE)
    message("==== Merging OTU and Sequence Tables ====")
    # Find path to data
    pathToSeqTables <- file.path(options$outDir, paste0(options$projectPrefix, "_sequence_tables"))
    seq.paths <- list.files(path = pathToSeqTables, recursive = TRUE,
                            pattern = "\\.rds$", 
                            full.names = TRUE)	# Get RDS files
    pathToOTUTables <- file.path(options$outDir, paste0(options$projectPrefix, "_taxonomy_tables"))
    otu.paths <- list.files(path = pathToOTUTables, recursive = TRUE,
                            pattern = "\\.rds$",
                            full.names = TRUE)  # Get CSV files
   
    #pathToOTUTables <- vector(mode="character", length=length(amplicons))
    #pathToSeqTables <- vector(mode="character", length=length(amplicons))
    track <- matrix(, nrow = length(amplicons), ncol = 0)
    #for (i in 1:length(amplicons)) {
    #  currTab <- amplicons[[i]]
    #  currSeq <- currTab[1]
    #  currOTU <- currTab[2]
    #  pathToOTUTables[i] <- currOTU
    #  pathToSeqTables[i] <- currSeq
    #}
    
    # Create final merged table
    #finalTable <- mergeSamples(unique(pathToOTUTables), unique(pathToSeqTables), options$projectPrefix, options$assignTaxLevels)
    finalTable <- mergeSamples(unique(otu.paths), unique(seq.paths), options$projectPrefix, options$assignTaxLevels)
    if (typeof(finalTable) == FALSE) {
      write("Unable to create merged table", log.file, append = TRUE)
      print("Unable to create final table")
    } else {
      write(paste0("Created: ", finalTable), log.file, append = TRUE)
    }
  }
   
  # Switch working directory back to original working directory
  setwd(curr.wd)
  
  print("Done! All outputs were sucessfully generated.")
  
}

# ======
# Analyze Single-End
# ======
help_single_analysis <- function(filtFs, sName, index, options, logFile) {
  # Dereplicate Sequences
  derepFs <- dada2::derepFastq(filtFs)
  
  # Learn Errors
  errFs <- dada2::learnErrors(filtFs, multithread = TRUE)
  
  # Denoise Data
  if (!is.null(options$dadaBandSize) & !is.null(options$dadaBandSize)) {
    dadaFs <- dada2::dada(derep = derepFs, err=errFs)
  } else {
    dadaFs <- dada2::dada(derepFs, err=errFs, BAND_SIZE=options$dadaBandSize, OMEGA_A=options$dadaOmegaA)
  }
  
  # Create table containing all unique sequences and their frequency within the original dataset
  seqtab <- dada2::makeSequenceTable(dadaFs, orderBy = "abundance")
  
  # Set .rds file name
  seq.print <- paste0(index, "_", sName, "_seqtab.rds")
  
  # Remove chimeras
  seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
  if (options$createChimeraDetectionTable) {
    seqtab.chim <- isBimeraDenovoTable(
      seqtab,
      minSampleFraction = options$chimeraDetectionMinSampleFraction,
      ignoreNNegatives = options$chimeraDetectionIgnoreNegatives,
      minFoldParentOverAbundance = options$chimeraDetectionMinFoldParentOverabundance,
      minParentAbundance = options$chimeraDetectionParentAbundance,
      allowOneOff = options$chimeraDetectionAllowOneOff,
      minOneOffParentDistance = options$chimeraDetectionMaxShift,
      maxShift = options$chimeraDetectionMaxShift,
      multithread = options$multithread,
      verbose = options$verbose
    )
    
    # Save Sequence Table with chimeras labeled
    chim.print <- saveSeqs(seqtab.chim, sName, index, options$outDir, add.table = TRUE, options$projectPrefix)
    write(paste0("Created: ", basename(chim.print)), logFile, append = TRUE)
  }
  
  # Save Sequence Table without chimeras
  seq.print <- saveSeqs(seqtab.nochim, sName, index, options$outDir, add.table = FALSE, options$projectPrefix)
  write(paste0("Created: ", basename(seq.print)), logFile, append = TRUE)
  
  # Classify sequences
  taxa <- dada2::assignTaxonomy(seqtab.nochim, options$taxDatabase, minBoot = options$assignTaxMinBootstrap, multithread=options$multithread)
  f.print <- saveTaxonomyTables(taxa, sName, options$outDir, index, options$projectPrefix)
  write(paste0("Created: ", basename(f.print)), logFile, append = TRUE)
  
  # Get counts
  derep <- sum(dada2::getUniques(derepFs))
  dada <- sum(dada2::getUniques(dadaFs))
  nochim <- rowSums(seqtab.nochim)
  
  # Combine all information
  track <- cbind(sName, derep, dada, nochim)
  
  # Return all required information
  return(c(seq.print, f.print, track))
  
}

single_analysis <- function(fp, sample.names, file, getQuality = FALSE, REGEX = "*.fastq(.gz)?$", logFile) {
  # Create object - simplifies debugging process
  object <- readConfig(file, isPaired = FALSE, type = c('auto', 'filter'))
  
  # Sort path to extract all .fastq files (external variables can be accessed)
  Fs <- sort(list.files(fp, pattern = REGEX, full.names = TRUE))
  
  # Filter and Trim (generates "filtered_objects.RData" file)
  write("==== Filtering and Trimming Single-End Amplicon Sequences ====", logFile, append = TRUE)
  if (!is.null(object)) {
    filtered.files <- filtTrim(sample.names=sample.names, object=object, forwardFs=Fs)
  } else {
    stop("Error created when creating filtering object")
  }
  
  # Plot Quality Distribution and save object
  if (getQuality) {
    write("==== Plotting quality distribution AFTER trimming ====", logFile, append = TRUE)
    plot_object <- readConfig(file, type = "qualityplot")
    if (!is.null(plot_object)) {
      plots <- plotQuality(filtered.files, object@projectPrefix, plot_object)
    } else {
      stop("Unable to plot object")
      
      # Write PDF files in output directory
      if (!is.null(plots)) {
        write(paste0("Created: ", paste0(object@projectPrefix, "_fastq2otu_quality_plots_AFTER.pdf")), logFile, append = TRUE)
        pdf(file = paste0(object@projectPrefix, "_fastq2otu_quality_plots_AFTER.pdf"))
        	print(plots)
        dev.off()
      } else {
        stop("Unable to write quality plots to PDF") # Should rarely execute. Mostly for debugging purposes
      }
    }
  }
   
  options <- yaml::yaml.load_file(file)
    
    
    # Finish analysis
    infoList <- lapply(1:length(Fs), function(input, samples, config) {
      return(help_single_analysis(filtered.files[[input]], sName=samples[input], index=input, options=config, logFile=logFile))
    }, samples=sample.names, config=options)
    
    message("\n==== Create Summary Table ====")
    sample_ls <- vector(mode="character", length=length(infoList))
    derep_ls <- vector(mode="character", length=length(infoList))
    dada_ls <- vector(mode="character", length=length(infoList))
    nochim_ls <- vector(mode="character", length=length(infoList))
    for (i in 1:length(infoList)) {
      currTab <- infoList[[i]]
      sample_ls[i] <- currTab[3]
      derep_ls[i] <- currTab[4]
      dada_ls[i] <- currTab[5]
      nochim_ls[i] <- currTab[6]
    }
    
    # Save summary table
    if (file.exists("filtered_objects.RData")) {
      load("filtered_objects.RData") # Load saveFilt obect in current environment
      
      # Combine all information
      track <- cbind(sample_ls, saveFilt, derep_ls, dada_ls, nochim_ls)
      
      # Show tracking and save to file
      s.print <- paste0(object@projectPrefix, "_read_retention_table.txt")
      write.table(track, file = s.print, sep = "\t")
      message("Created: ", s.print)
    } else {
      stop("File Not Found Error: Could not find filtered_objects.RData file")
    }
    
    # Return final information
    return(infoList)
    
  }

# ======
# Analyze Paired-End
# ======
help_paired_analysis <- function(filtFs, filtRs, sName, index, options) {
    if (is.null(filtFs)) {
      stop("Invalid input for filtFs.")
    }
    
    # Dereplicate Sequences
    derepFs <- dada2::derepFastq(filtFs)
    derepRs <- dada2::derepFastq(filtRs)
    
    # Learn Errors
    errFs <- dada2::learnErrors(filtFs, multithread = TRUE)
    errRs <- dada2::learnErrors(filtRs, multithread = TRUE)
    
    # Denoise Data
    if (!is.null(options$dadaOmegaA) & !is.null(options$dadaBandSize)) {
      dadaFs <- dada2::dada(derep = derepFs, err=errFs)
      dadaRs <- dada2::dada(derep = derepRs, err=errRs)
    } else {
      dadaFs <- dada2::dada(derepFs, err=errFs, BAND_SIZE=options$dadaBandSize, OMEGA_A=options$dadaOmegaA)
      dadaRs <- dada2::dada(derepRs, err=errRs, BAND_SIZE=options$dadaBandSize, OMEGA_A=options$dadaOmegaA)
    }
    
    # Merge pairs
    merged_amplicons <- mergeSeqPairs(dadaFS=dadaFs, dadaRS=dadaRs, derepFS=derepFs, derepRS=derepRs, options)
    if (is.null(merged_amplicons)) {
      message("Unable to merge sequences\n")
      return(c(errFs, errRs)) # Return error objects without proceeding forward
    }
    
    # Create table containing all unique sequences and their frequency within the original dataset
    seqtab <- dada2::makeSequenceTable(merged_amplicons, orderBy = "abundance")
    
    seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
    if (options$createChimeraDetectionTable) {
      seqtab.chim <- isBimeraDenovoTable(
        seqtab,
        minSampleFraction = options$chimeraDetectionMinSampleFraction,
        ignoreNNegatives = options$chimeraDetectionIgnoreNegatives,
        minFoldParentOverAbundance = options$chimeraDetectionMinFoldParentOverabundance,
        minParentAbundance = options$chimeraDetectionParentAbundance,
        allowOneOff = options$chimeraDetectionAllowOneOff,
        minOneOffParentDistance = options$chimeraDetectionMaxShift,
        maxShift = options$chimeraDetectionMaxShift,
        multithread = options$multithread,
        verbose = options$verbose
      )
      
      # Save Sequence Table with chimeras labeled
      chim.print <- saveSeqs(seqtab.chim, sName, index, options$outDir, add.table = TRUE, options$projectPrefix)
    }
    
    # Save Sequence Table without chimeras
    seq.print <- saveSeqs(seqtab.nochim, sName, index, options$outDir, add.table = FALSE, options$projectPrefix)
    
    # Classify sequences
    taxa <- dada2::assignTaxonomy(seqtab.nochim, options$taxDatabase, minBoot = options$assignTaxMinBootstrap, multithread=options$multithread)
    
    # Save OTU table
    f.print <- saveTaxonomyTables(taxa, sName, options$outDir, index, options$projectPrefix)
    
    # Get counts
    forward.derep <- sum(dada2::getUniques(derepFs))
    reverse.derep <- sum(dada2::getUniques(derepRs))
    forward.dada <- sum(dada2::getUniques(dadaFs))
    reverse.dada <- sum(dada2::getUniques(dadaRs))
    merged <- sum(dada2::getUniques(merged_amplicons))
    nochim <- rowSums(seqtab.nochim)
    
    # Combine all information
    track <- cbind(sName, forward.derep, reverse.derep, forward.dada, reverse.dada, merged, nochim)
    
    # Return all required information
    return(c(seq.print, f.print, track))
  }
  
paired_analysis <- function(fp, sample.names, file, getQuality = FALSE, REGEX_1 = "*_1.fastq(.gz)?$", REGEX_2 = "*_2.fastq(.gz)?$", logFile) {
    # Sort path to extract all .fastq files (allows .fastq or .fastq.gz file extensions)
    Fs <- sort(list.files(fp, pattern = REGEX_1, full.names = TRUE))
    Rs <- sort(list.files(fp, pattern = REGEX_2, full.names = TRUE))
    
    if (length(Fs) != length(Rs)) {
      stop("Error: Unequal number of files containing forward and reverse sequences")
    }
    
    # Create filter and trim object
    object <- readConfig(file, isPaired = TRUE, type = c('auto', 'filter'))
    
    # Filter and Trim
    write("==== Filtering and Trimming Paired-End Amplicon Sequences ====", logFile, append = TRUE)
    if (!is.null(object)) {
      filtered.files <- filtTrim(sample.names=sample.names, object=object, forwardFs=Fs, reverseRs=Rs)
    } else {
      stop("Error generated when creating fastFilt object")
    }
    
    # Create output file names for filtered sequences
    filt.forward <- filtered.files[grep("*_R1_filt_trimmed.fastq.gz", filtered.files)]
    filt.reverse <- filtered.files[grep("*_R2_filt_trimmed.fastq.gz", filtered.files)]
    
    # Plot Quality Distribution and save object
    if (getQuality) {
      write("==== Plotting quality distribution AFTER trimming ====", logFile, append = TRUE)
      plot_object <- readConfig(file, type = "qualityplot")
      ordered.files <- c(rbind(filt.forward, filt.reverse)) # Alternate plots from forward and reverse sequences
      plots <- plotQuality(ordered.files, options$projectPrefix, plot_object)
      
      # Write PDF files in output directory
      if (!is.null(plots)) {
        write(paste0("Created: ", paste0(object@projectPrefix, "_fastq2otu_quality_plots_AFTER.pdf")), logFile, append = TRUE)
        pdf(file = paste0(object@projectPrefix, "_fastq2otu_quality_plots_AFTER.pdf"))
        	print(plots)
        dev.off()
      } else {
        stop("Unable to write quality plots to PDF") # Should rarely execute. Mostly for debugging purposes
      }
    }
    
    # Create object - simplifies debugging process -- TODO: Test all inputs in config file
    options <- yaml::yaml.load_file(file)
    infoList <- lapply(1:length(Fs), function(input, samples, label, config) {
      return(help_paired_analysis(filt.forward[[input]], filt.reverse[[input]], sName=samples[input], index=input, options=config))
    },  samples=sample.names, label=1:length(sample.names), config=options)
    
    # Get summary table
    message("\n==== Create Summary Table ====")
    sample_ls <- vector(mode="character", length=length(infoList))
    derepFs_ls <- vector(mode="character", length=length(infoList))
    derepRs_ls <- vector(mode="character", length=length(infoList))
    dadaFs_ls <- vector(mode="character", length=length(infoList))
    dadaRs_ls <- vector(mode="character", length=length(infoList))
    merged_ls <- vector(mode="character", length=length(infoList))
    nochim_ls <- vector(mode="character", length=length(infoList))
    for (i in 1:length(infoList)) {
      currTab <- infoList[[i]]
      sample_ls[i] <- currTab[3]
      derepFs_ls[i] <- currTab[4]
      derepRs_ls[i] <- currTab[5]
      dadaFs_ls[i] <- currTab[6]
      dadaRs_ls[i] <- currTab[7]
      nochim_ls[i] <- currTab[9]
      merged_ls[i] <- currTab[8]
    }
    
    # Save summary table
    if (file.exists("filtered_objects.RData")) {
      load("filtered_objects.RData") # Load saveFilt obect in current environment
      
      # Combine all information
      track <- cbind(sample_ls, saveFilt, derepFs_ls, derepRs_ls, dadaFs_ls, dadaRs_ls, nochim_ls, merged_ls)
      
      # Show tracking and save to file
      s.print <- paste0(object@projectPrefix, "_read_retention_table.txt")
      write.table(track, file = s.print, sep = "\t")
      message("Created: ", s.print)
    } else {
      stop("File Not Found Error: Could not find filtered_objects.RData file")
    }
    
    # Return final information
    return(infoList)
  }
