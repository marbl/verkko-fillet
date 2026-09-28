#!/bin/bash

verkko_original_dir=$(realpath $1)
verkko_fillet_dir=$(realpath $2)

if [ -d "$verkko_fillet_dir" ]; then
    rm -rf "$verkko_fillet_dir"
fi

for folder in 5-untip 6-rukki 8-hicPipeline
do
if [ -d $verkko_original_dir/$folder ]; then
mkdir -p $verkko_fillet_dir/$folder
fi
done

# Change directory to the target directory
cd "$verkko_fillet_dir" || { echo "Error: Cannot change directory to $verkko_fillet_dir"; exit 1; }

echo -e "verkko_original_dir : $verkko_original_dir"
echo -e "verkko_fillet_dir : $verkko_fillet_dir"

# Array of files and directories to link
files=("assembly.fasta" "assembly.scfmap" "assembly.colors.csv" "assembly.homopolymer-compressed.gfa" "assembly.homopolymer-compressed.noseq.gfa" "assembly.paths.tsv" "assembly.fasta.fai" "8-hicPipeline/rukki.paths.gaf" "5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.hifi-coverage.csv" "8-hicPipeline/unitigs.hpc.fasta" "5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.gfa" "8-hicPipeline/hic.byread.compressed" "8-hicPipeline/final_contigs/assembly.ont-coverage.csv" "8-hicPipeline/final_contigs/assembly.hifi-coverage.csv")

# Iterate over the files and create symbolic links, handling subdirectories recursively
for file in "${files[@]}"; do
    src="$verkko_original_dir/$file"
    dest="$verkko_fillet_dir/$file"
    if [ -e "$src" ]; then
        # Create parent directory if it doesn't exist
        mkdir -p "$(dirname "$dest")"
        cmd="ln -s \"$src\" \"$dest\""
        echo $cmd && eval $cmd
    else
        echo "Error: $file does not exist in $verkko_original_dir."
    fi
done


# ln fasta
folders=("6-rukki" "5-untip")  # Array of folder names

for folder in "${folders[@]}"; do
    echo "Processing folder: $folder"

    if ls "$verkko_original_dir/$folder/unitig"*.fasta >/dev/null 2>&1; then
        # If files exist, create symbolic links for each
        for file in "$verkko_original_dir/$folder/unitig"*.fasta; do
            cmd="ln -s \"$file\" \"$folder/\""
            echo $cmd && eval $cmd
        done
    else
        echo "Error: No files matching the pattern 'unitig*.fasta' exist in $verkko_original_dir/$folder."
    fi
done
