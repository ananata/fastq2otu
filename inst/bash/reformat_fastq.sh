#!/bin/bash
################################################################################
# reformat_fastq.sh
#
# Seperates interleaved fastq file into two seperate files. Please verify that
# file is indeed interleaved before executing the script.
#
# usage
#   reformat_fastq.sh -i ../raw/*.gz -o ../path/to/output
#
# Adapted from https://gist.github.com/nathanhaigh/3521724
# Inspired by Torsten Seemann's blog post:
# http://thegenomefactory.blogspot.com.au/2012/05/cool-use-of-unix-paste-with-ngs.html
################################################################################

usage(){
echo "
Modified by Nana Afia Twumasi-Ankrah for FASTQ2OTU R package
Last Updated: August 11, 2020
Adapted from https://gist.github.com/nathanhaigh/3521724
Inspired by Torsten Seemann's blog post: http://thegenomefactory.blogspot.com.au/2012/05/cool-use-of-unix-paste-with-ngs.html

Input parameters:
-i | --input =<file>         Main input. Interleaved fastq file
-o | --output=<dir>          Path to output directory.

}

while test $# -gt 0; do
           case "$1" in
                -i|--input)
                    shift
                    first_argument=$1
                    shift
                    ;;
                -o|--output)
                    shift
                    last_argument=$1
                    shift
                    ;;
                *)
                   echo "$1 is not a recognized flag!"
                   return 1;
                   ;;
          esac
  done

  echo "Input File : $first_argument";
  echo "Output Directory : $last_argument";


filename=$(basename -- "$first_argument")
extension="${filename##*.}"
base="${filename%.*}"
suffix=".fastq"

echo "Basename: $base"

for f in $first_argument; do
    echo "Splitting $f"

    r1=${last_argument}/${base}"_1.fastq"
    r2=${last_argument}/${base}"_2.fastq"
    echo "Creating $r1..."
    echo "Creating $r2..."

    # Split the file into two
    cat $f | paste - - - - - - - - | tee >( \
    cut -f 1-4 | \tr "\t" "\n" | pigz --best >$r1 \
    ) | \
    cut -f 5-8 | tr "\t" "\n" | pigz --best > $r2
done
