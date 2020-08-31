#!/usr/bin/env bash
# A shell script to download fastq files from NCBI's SRA database
# Written by: Nana Afia Twumasi-Ankrah
# Last updated on: 08/30/2020
# -------------------------------------------------------

# ERROR CHECKS:
# Verify the type of input and number of values
# Display an error message if no inputs are not provided
# Exit the shell script with a status of 1 using exit 1 command.

# Check to see if any arguments are provided
[ $# -eq 0 ] && { echo "Usage: $0 path/to/fastq-dump path/to/SRR_Acc_List.txt path/to/output/"; exit 1; }


fastq=$1 # Path to fastq-dump script (keep in same directory)
[[ ! -f $fastq ]] && { echo "fastq-dump script could not be found (invalid path): $fastq"; exit 1; }

input=$2 # Path to list of SRA ids
[[ ! -f $input ]] && { echo "Input file could not found (invalid path): $input"; exit 1; }

output=$3 # Path to output directory


# Create output directory if it doesn't already exist
if [ ! -d "$output" ]; then
	read -p "Are you sure you want to create a new directory? <y/N> ";
	if [[ $prompt == "y" || $prompt == "Y" || $prompt == "yes" || $prompt == "Yes" ]]
	then
 		echo "Creating new directory: $output";
		mkdir $output;
	else
		echo "Unable to create new directory. Please enter new path.";
  		exit 0;
	fi
fi

# Navigate to output directory
cd $output

# Read list and install files
while IFS= read -r line
do
	$fastq --split-3 --gzip "$line"
done < "$input"

