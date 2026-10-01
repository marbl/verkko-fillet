#!/bin/bash
# rDNA=/vf/users/Phillippy/projects/giraffeT2T/assembly/verkko2.2_hifi-duplex_trio-hic/rDNA_seq/NT_167214.1.fa
rDNA=$1
cpu=$2
fa=$3

# load if module system is available
if command -v module &> /dev/null
then
    module load mashmap
fi
# check if mashmap is installed
if ! command -v mashmap &> /dev/null
then
    echo "mashmap could not be found, please install mashmap first."
    exit 1
fi

mkdir -p rDNA_mapping
mashmap -t $cpu --noSplit -q $rDNA \
	-r $fa -s 13357 --pi 85 -f none -o rDNA_mapping/45S.mashmap.out

cat rDNA_mapping/45S.mashmap.out |\
awk -F'\t' '{
  split($13, id, ":");
  print $6, $8, $9, "45S", (id[3]*100), $5;
}' OFS='\t' > rDNA_mapping/45S.mashmap.bed