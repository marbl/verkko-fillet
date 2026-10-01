#!/bin/bash
# Run with bash

# Require a fasta path as the first argument
if [ -z $1 ]; then
	echo "Usage: ./find_telomere.sh <fasta>"
	exit -1
fi

# Store the input fasta path passed as the first argument
file=$1
# Strip the directory portion, keeping only the filename
file_name=$(basename $file)
# VGP_PIPELINE=$(readlink -f "$0")
# Resolve the vgp-assembly toolkit directory relative to this script's location
VGP_PIPELINE=$(realpath $(dirname "$0"))/vgp-assembly/

# Symlink the fasta into the current directory if it isn't already here
if [ ! -e $file_name ]; then
	ln -s $file
fi

# Work with the local (symlinked) copy from here on
file=$file_name
# Derive an output prefix by dropping a .fasta or .fa extension
prefix=$(echo $file | sed 's/.fasta$//g' | sed 's/.fa$//g')

# Call the telomere finder and reformat its output to contig/coords/score columns
$VGP_PIPELINE/telomere/find_telomere $file | awk '{print $1"\t"$(NF-4)"\t"$(NF-3)"\t"$(NF-2)"\t"$(NF-1)"\t"$NF}' - > $prefix.telomere
# Mask low-complexity regions in the fasta
sdust $file > $prefix.sdust
# Compute per-sequence lengths from the fasta
java -cp $VGP_PIPELINE/telomere/telomere.jar SizeFasta $file > $prefix.lens

# Lowering threshold to 0.10 (10%) from the initial 0.40 (40%)
# Slide windows over the telomere calls to flag telomere-rich regions
java -cp $VGP_PIPELINE/telomere/telomere.jar FindTelomereWindows $prefix.telomere 99.9 0.1 > $prefix.windows
# Identify candidate assembly breaks using lengths, low-complexity, and telomere calls
java -cp $VGP_PIPELINE/telomere/telomere.jar FindTelomereBreaks $prefix.lens $prefix.sdust $prefix.telomere > $prefix.breaks