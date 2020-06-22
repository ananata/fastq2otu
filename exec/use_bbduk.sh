#!/bin/bash
# Script location: /dada2/resources/use_bbduk.sh
# Usage: ./use_bbduk.sh /path/to/raw_fastq /path/to/noprimer_seqs /path/to/sample_ids.txt /path/to/adapters.fa
# Author: Nana Afia Twumasi-Ankrah
# Date Posted: 03/16/2020
#
# Use BBDuk tool (apart of BBTools package) to remove adapter sequences
# Get more information about the tool at: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/

# Keep in same directory
bbduk="bbduk.sh"

# Set variables
INDIR= $1 # Name and path of directory containing raw fastq files
OUTDIR= $2 # Name and path of directory to create
ADAPTERS=$4 # Set path to file containing adapters
ISPAIRED=$5 # If value is present, then the data is paired, if no value present the data is single-end
# Create output directory
if [ ! -d $OUTDIR ]
then
	mkdir $OUTDIR
fi

# Get list of samples
SAMPLES=$3

# Record parameters
echo "Current date: $(date)" >> ${OUTDIR}/bbduk_parameters.txt
echo ${INDIR} >> ${OUTDIR}/bbduk_parameters.txt
echo ${OUTDIR} >> ${OUTDIR}/bbduk_parameters.txt
echo ${ADAPTERS} >> ${OUTDIR}/bbduk_parameters.txt

# Change name to correspond to current project and move to project folder
# in1 = Path to forward read(s)
# in2 = Path to reverse read(s)
# out1/2 = Path(s) to output files
# ref = Path to FASTA-formatted file containing adapter sequences
# ktrim = Read/trim kmers from the left (beginning of sequence)
# k = Length of adapters (scans sequences K-nucleotides at a time)
# tpo = Trims adapters based on pair overlap detection using BBMerge
# tpe = Trims both reads to the same length (in the event that an adapter kmer was only detected in one of them).
# mink = Minimum k-value at read end
# hdist = Number of mismatches allowed when searching query against reference.
if [ ! -z "$5" ] # Checks to see if fifth argument is present
then
	for sample in $(cat $SAMPLES)
	do
		echo "On sample: $sample"
        	# Assumes equal formatting across files (feel free to modify to suit file format)
        	${bbduk} in=${INDIR}/${sample}.fastq out=${OUTDIR}/${sample}_no_adapter.fq.gz ref=${ADAPTERS} mink=11 hdist=1 ktrim=l k=33 tbo tpe >> ${OUTDIR}/${sample}_adapter_trimming_summary.txt 2>&1
	done
else
	for sample in $(cat $SAMPLES)
	do
        	echo "On sample: $sample"
        	# Assumes equal formatting across files (feel free to modify to suit file format)
		${bbduk} in1=${INDIR}/${sample}_1.fastq in2=${INDIR}/${sample}_2.fastq out1=${OUTDIR}/${sample}_R1_no_adapter.fq.gz out2=${OUTDIR}/${sample}_R2_no_adapter.fq.gz ref=${ADAPTERS} mink=11 hdist=1 ktrim=l k=33 tbo tpe >> ${OUTDIR}/${sample}_adapter_trimming_summary.txt 2>&1
	done
fi
