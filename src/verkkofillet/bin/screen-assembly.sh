#!/bin/bash

print_help() {
        cat <<'EOF'
Usage:
    screen-assembly.sh <verkko_dir> [threads] [minlength] [ebv_fasta] [rdna_fasta] [mt_fasta]

Arguments:
    verkko_dir   Path to a Verkko output directory.
    threads      Number of threads (default: 10).
    minlength    Minimum contig length for screening (default: 1000).
    ebv_fasta    EBV contaminant FASTA (default: Verkko bundled EBV FASTA).
    rdna_fasta   rDNA contaminant FASTA (default: Verkko bundled rDNA FASTA).
    mt_fasta     MT contaminant FASTA (default: Verkko bundled MT FASTA).

Help:
    screen-assembly.sh helpme
    screen-assembly.sh -h
    screen-assembly.sh --help
EOF
}

if [ "$1" = "helpme" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
        print_help
        exit 0
fi

if [ -z "$1" ]; then
        print_help
        exit 1
fi

# check if module available and load samtools and seqtk
if command -v module &> /dev/null
then
    module load samtools
    module load seqtk
    ml perl
    ml verkko
fi
# check seqtk and mashamp is installed
tools="seqtk mashmap perl verkko"
for tool in $tools; do
    if ! command -v $tool &> /dev/null
    then
        echo "$tool could not be found, please install $tool first."
        exit 1
    fi
done

verkko=$(dirname $(which verkko))
verkko_script=$verkko/../lib/verkko/scripts/
verkko_data=$verkko/../lib/verkko/data

verkko_dir=$1
outPrefix=assembly_screen

if [ ! -d "$verkko_dir" ]; then
    echo "Verkko directory $verkko_dir does not exist."
    exit 1
fi

if [ -n "$2" ]; then
    threads=$2
else
    threads=10
fi

# if $3 is provided, use it as the minlength, otherwise default to 1000
if [ "$3" == "None" ]; then
    minlength=1000
else
    minlength=$3
fi

if [ -n "$4" ] && [ "$4" != "None" ]; then
    ebv=$4
else
    ebv=$verkko_data/human-ebv-AJ507799.2.fasta.gz
fi

if [ -n "$5" ] && [ "$5" != "None" ]; then
    rDNA=$5
else
    rDNA=$verkko_data/human-rdna-KY962518.1.fasta.gz
fi

if [ -n "$6" ] && [ "$6" != "None" ]; then
    mt=$6
else
    mt=$verkko_data/human-mito-NC_012920.1.fasta.gz
fi


echo -e "Screening assembly in $verkko_dir for contaminants:\nEBV: $ebv\nrDNA: $rDNA\nMT: $mt\nMinimum contig length: $minlength\nThreads: $threads"

for i in $verkko_script/screen-assembly.pl $verkko_dir/assembly.fasta $verkko_dir/assembly.homopolymer-compressed.noseq.gfa $verkko_dir/assembly.scfmap $verkko_dir/5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.hifi-coverage.csv $ebv $rDNA $mt; do
    if [ ! -f "$i" ]; then
        echo "Error: Required file $i does not exist."
        exit 1
    fi
done


cmd="perl $verkko_script/screen-assembly.pl \
--assembly $verkko_dir/assembly.fasta \
--threads $threads \
--graph $verkko_dir/assembly.homopolymer-compressed.noseq.gfa \
--graphmap $verkko_dir/assembly.scfmap \
--hifi-coverage $verkko_dir/5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.hifi-coverage.csv \
--minlength $minlength \
--output  $outPrefix \
--contaminant EBV $ebv \
--contaminant rDNA $rDNA \
--contaminant MT $mt"

echo -e $cmd
eval $cmd