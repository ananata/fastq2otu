#!/bin/bash
# A shell script to download fastq files from NCBI's SRA database
# Written by: Nana Afia Twumasi-Ankrah
# Last updated on: 03/20/2020
# -------------------------------------------------------

# Verify the type of input and number of values
# Display an error message if both inputs are not provided
# Exit the shell script with a status of 1 using exit 1 command.

[ $# -eq 0 ] && { echo "Usage: $0 path/to/fastq-dump path/to/SRR_Acc_List.txt path/to/output/"; exit 1; }


fastq="fastq-dump" # Path to fastq-dump script (keep in same directory)
input=$2 # Path to list of SRA ids
output=$3 # Path to output directory


# Create output directory
mkdir $output
cd $output

# Read list and install files
while IFS= read -r line
do
  $fastq --split-3 --gzip "$line"
done < "$input"

