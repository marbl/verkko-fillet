#!/bin/bash

rDNA_gap=$1
fasta=$2
out_fasta=$3
gap_size=$4
force=$5 #True or False

log=rDNA_gap_cleaning.log

if [ "$force" = "True" ]; then
    echo "Force flag is set. Overwriting existing output file."
    rm -rf $out_fasta.gz
fi
if [ "$force" != "True" ] && [ -f "$out_fasta.gz" ]; then
    echo "Error : Output fasta file $out_fasta.gz already exists. Use force flag to overwrite."
    exit 1
fi


# check if module available and load samtools and seqtk 
if command -v module &> /dev/null
then
    module load samtools
    module load seqtk
fi

# check if samtools and seqtk is installed
if ! command -v samtools &> /dev/null
then
    echo "samtools could not be found, please install samtools first."
    exit
fi

if ! command -v seqtk &> /dev/null
then
    echo "seqtk could not be found, please install seqtk first."
    exit
fi

# check input files exist
if [ ! -f "$rDNA_gap" ]; then
    echo "rDNA gap file $rDNA_gap does not exist."
    exit 1
fi
if [ ! -f "$fasta" ]; then
    echo "Fasta file $fasta does not exist."
    exit 1
fi

if [ ! -f "$fasta.fai" ]; then
    echo "Fasta index file $fasta.fai does not exist. Creating index..."
    samtools faidx $fasta
fi


(echo ">rDNA_gap"; head -c $gap_size /dev/zero | tr '\0' 'N') | seqtk seq -l 60 - > rDNA_gap.fasta


print_progress() {
    current=$1
    total=$2
    width=40

    if [ "$total" -le 0 ]; then
        return
    fi

    percent=$(( current * 100 / total ))
    filled=$(( current * width / total ))
    empty=$(( width - filled ))

    bar=$(printf "%${filled}s" "" | tr ' ' '#')
    space=$(printf "%${empty}s" "")
    printf "\rProgress: [${bar}${space}] %3d%% (%d/%d)" "$percent" "$current" "$total"

    if [ "$current" -eq "$total" ]; then
        printf "\n"
    fi
}

if [ -f "$out_fasta" ]; then
    echo "Output fasta file $out_fasta already exists. Overwriting."
    rm -f "$out_fasta"
fi

contig_list_tmp=$(mktemp)
cut -f 1 $fasta.fai > "$contig_list_tmp"
total_contigs=$(wc -l < "$contig_list_tmp")
processed=0
: > "$out_fasta"

while IFS= read -r contig; do
    [ -z "$contig" ] && continue
    echo "Processing contig: $contig"
    isingap=$(awk -F '\t' -v c="$contig" 'NR > 1 && $1 == c && $4 == "True" {n++} END {print n+0}' "$rDNA_gap")
    if [ $isingap -eq 0 ]; then
        # echo "No gaps found for contig $contig in $rDNA_gap. Skipping."
        samtools faidx "$fasta" "$contig" >> "$out_fasta"
    elif [ $isingap -eq 1 ]; then
        echo "One gap found for contig $contig in $rDNA_gap. Extracting sequence."
        gap_start=$(awk -F '\t' -v c="$contig" 'NR > 1 && $1 == c && $4 == "True" {print $2; exit}' "$rDNA_gap")
        gap_end=$(awk -F '\t' -v c="$contig" 'NR > 1 && $1 == c && $4 == "True" {print $3; exit}' "$rDNA_gap")
        if [ -z "$gap_start" ] || [ -z "$gap_end" ]; then
            echo "Invalid gap coordinates for $contig. Writing original sequence."
            samtools faidx "$fasta" "$contig" >> "$out_fasta"
            processed=$((processed + 1))
            print_progress "$processed" "$total_contigs"
            continue
        fi
        echo "Gap coordinates: $gap_start-$gap_end"
        # Extract the sequence from the fasta file using seqtk
        samtools faidx "$fasta" "$contig:1-$gap_start" > tmp1.fasta || { echo "Failed to extract left segment for $contig"; rm -f tmp1.fasta tmp2.fasta; continue; }
        samtools faidx "$fasta" "$contig:$gap_end-" > tmp2.fasta || { echo "Failed to extract right segment for $contig"; rm -f tmp1.fasta tmp2.fasta; continue; }
        (echo ">$contig"; cat tmp1.fasta rDNA_gap.fasta tmp2.fasta | \
        grep -v '^>' | tr -d '\n') | \
        seqtk seq -l 60 - >> "$out_fasta"
        rm -f tmp1.fasta tmp2.fasta
        echo -e -n "$contig has a gap. Replaced with rDNA_gap sequence.\n" >> "$log"
    else
        echo "Multiple gaps found for contig $contig in $rDNA_gap. Skipping."
        gap_starts=$(awk -F '\t' -v c="$contig" 'NR > 1 && $1 == c && $4 == "True" {print $2}' "$rDNA_gap" | sort -n | head -1)
        gap_ends=$(awk -F '\t' -v c="$contig" 'NR > 1 && $1 == c && $4 == "True" {print $3}' "$rDNA_gap" | sort -nr | head -1)
        if [ -z "$gap_starts" ] || [ -z "$gap_ends" ]; then
            echo "Invalid multi-gap coordinates for $contig. Writing original sequence."
            samtools faidx "$fasta" "$contig" >> "$out_fasta"
            processed=$((processed + 1))
            print_progress "$processed" "$total_contigs"
            continue
        fi
        echo "Using first gap start: $gap_starts and last gap end: $gap_ends for contig $contig."
        samtools faidx "$fasta" "$contig:1-$gap_starts" > tmp1.fasta || { echo "Failed to extract left segment for $contig"; rm -f tmp1.fasta tmp2.fasta; continue; }
        samtools faidx "$fasta" "$contig:$gap_ends-" > tmp2.fasta || { echo "Failed to extract right segment for $contig"; rm -f tmp1.fasta tmp2.fasta; continue; }
        (echo ">$contig"; cat tmp1.fasta rDNA_gap.fasta tmp2.fasta | \
        grep -v '^>' | tr -d '\n') | \
        seqtk seq -l 60 - >> "$out_fasta"
        rm -f tmp1.fasta tmp2.fasta
        echo -e "$contig has multiple gaps. Used first gap start: $gap_starts and last gap end: $gap_ends for replacement." >> "$log"
    fi
    processed=$((processed + 1))
    print_progress "$processed" "$total_contigs"
done < "$contig_list_tmp"

rm -f "$contig_list_tmp"

rm -f rDNA_gap.fasta

echo "Completed processing all contigs. Output written to $out_fasta."


bgzip "$out_fasta" && echo "Compressed $out_fasta to $out_fasta.gz"

echo -e "Indexing $out_fasta.gz..."
samtools faidx "$out_fasta.gz"