#!/bin/bash

# this script will take a reference (or assume human) and map the graph as well as final assembly to it
# assumes HPC version of reference is named same but .hpc.fasta, if it doesn't exist it will be created
# the output is a csv file suitable for bandage to new chromosome assignments
# it will also output assignments of the final sequences to each chromosome

[ -n "$DEBUG" ] && set -x

p=`pwd`

# make help
if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
   echo "Usage: $0 [reference.fasta] [identity_threshold] [contigs_to_chr_names.tsv] [chr_name_for_contigs]"
   echo "  reference.fasta: path to reference fasta file (default is human T2T CHM13 v1.1)"
   echo "  identity_threshold: minimum percent identity for a contig to be assigned to a chromosome (default is 99)"
   echo "  contigs_to_chr_names.tsv: optional tsv file mapping contig names to desired chromosome names. If not provided, contigs will be assigned to chromosomes based on best mashmap hit to reference. Format: <contig_name_in_fasta> <desired_chromosome_name>"
   echo "  chr_name_for_contigs: optional chromosome name to assign to contigs that do not have a strong hit to the reference (default is . for unassigned)"
   echo "  Examples:"
   echo "    $0 chr1.fasta 99 contigs_to_chr_names.tsv chr1"
   echo "    $0 chr1.fasta 99 . ."
   exit 0
fi

# check dependency
if [[ "$1" == "check_dependencies" ]]; then
   echo "Checking dependencies..."
   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

   missing=0
   for tool in mashmap seqtk; do
      if ! command -v "$tool" >/dev/null 2>&1; then
         echo "Error: $tool not found on \$PATH"
         missing=1
      fi
   done

   neigh=$(which neighborhood 2>/dev/null)
   if [ -z "$neigh" ] && [ -x "$SCRIPT_DIR/neighborhood" ]; then
      neigh="$SCRIPT_DIR/neighborhood"
   fi
   if [ -z "$neigh" ] || [ ! -x "$neigh" ]; then
      echo "Error: gfacpp 'neighborhood' not found (\$PATH or $SCRIPT_DIR/neighborhood)"
      missing=1
   fi

   if [ "$missing" -ne 0 ]; then
      exit 1
   fi
   echo "All dependencies found"
   exit 0
fi


mashmap=$(which mashmap 2>/dev/null)
if [ "x$mashmap" == "x" ]; then
   echo "Error: mashmap not found"
   exit 1
fi

seqtk=$(which seqtk 2>/dev/null)
if [ "x$seqtk" == "x" ]; then
   echo "Error: seqtk not found"
   exit
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

neigh=$(which neighborhood 2>/dev/null)
if [ "x$neigh" == "x" ] && [ -x "$SCRIPT_DIR/neighborhood" ]; then
   neigh="$SCRIPT_DIR/neighborhood"
fi
if [ ! -x "$neigh" ]; then
   echo "Error: gfacpp 'neighborhood' not found"
   exit 1
fi


echo "Using mashmap: $mashmap"
echo "Using seqtk: $seqtk"
echo "Using neighborhood: $neigh"



if [[ "$#" -ge 1 ]]; then
   echo "Using custom reference sequence $1"
   ref=$(realpath -- "$1")
fi

if [ -z "$ref" ]; then
   echo "Error: reference fasta is required (run with -h for usage)"
   exit 1
fi

hpcRef=`echo "$ref"|sed s/.fasta$/.hpc.fasta/g |sed s/.fa$/.hpc.fasta/g`

if [[ "$ref" == *"hpc.fa"* ]]; then
   echo "$1 already hpc-ed"
   hpcRef=$ref
   echo $hpcRef
else
   if [ ! -e $hpcRef ]; then
      $seqtk hpc $ref > $hpcRef
   fi
fi

idy=99
if [[ "$#" -ge 2 ]] ; then
   echo "Using custom identity $2"
   idy=$2
fi
contigs="assembly.fasta"
if [[ "$#" -ge 3 ]] ; then
   echo "Using custom contigs file $3"
   contigs=$3
fi
if [ $idy -lt 0 ] || [ $idy -gt 100 ] ; then
   echo "Error: identity must be between 0-100 and you provided $idy"
   exit 1
fi

if [ ! -e assembly.mashmap.out ]; then
   # $mashmap -r $ref -q $contigs --pi 95 -s 10000 -t 8  -o >(awk '$11 >= 50000' > assembly.mashmap.out)
   "$mashmap" -r "$ref" -q "$contigs" --pi 95 -s 10000 -t 8 \
        -o temp.mashmap.out
    awk '$11 >= 50000' temp.mashmap.out > assembly.mashmap.out
    rm temp.mashmap.out
fi


utigs_mashmap=unitigs.hpc.mashmap.out

# Detect which Verkko sub-pipeline produced the unitigs (used both for picking
# the mashmap input below and for choosing the noseq.gfa to symlink later).
isRUK=`ls 6-rukki/unitig*.fasta 2>/dev/null |wc -l |awk '{print $1}'`
isHIC=`ls 8-hicPipeline/unitigs.hpc.fasta 2>/dev/null |wc -l |awk '{print $1}'`

if [ -e contigs.gfa ]; then
   cat contigs.gfa |awk '{if (match($1, "^S")) print $1"\t"$2"\t*\tLN:i:"length($3); else print $0}'  > tmp.gfa
   /data/korens/devel/verkko-tip/lib/verkko/scripts/inject_coverage.py --allow-absent 5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.hifi-coverage.csv tmp.gfa > contigs.noseq.gfa
   rm tmp.gfa
fi

if [ ! -e $utigs_mashmap ]; then
   isRUK=`ls 6-rukki/unitig*.fasta 2>/dev/null |wc -l |awk '{print $1}'`
   isHIC=`ls 8-hicPipeline/unitigs.hpc.fasta 2>/dev/null |wc -l |awk '{print $1}'`

   if [ $isRUK -gt 0 ]; then
      cd 6-rukki
      # $mashmap -r $hpcRef -q unitig*.fasta --pi 95 -s 10000 -f none -t 8 -o >(awk '$11 >= 50000' > unitigs.hpc.mashmap.out) 
      if [ ! -s unitigs.hpc.mashmap.raw ]; then
         $mashmap -r "$hpcRef" -q unitig*.fasta --pi 95 -s 10000 -f none -t 8 -o unitigs.hpc.mashmap.raw
      fi
      # awk '$11 >= 10000' unitigs.hpc.mashmap.raw > unitigs.hpc.mashmap.out
      awk '$11 >= 50000' unitigs.hpc.mashmap.raw > unitigs.hpc.mashmap.out

      cd ..
      ln -s 6-rukki/unitigs.hpc.mashmap.out $utigs_mashmap
   elif [ $isHIC -gt 0 ]; then
      cd 8-hicPipeline
      if [ ! -s unitigs.hpc.mashmap.raw ]; then
         # $mashmap -r $hpcRef -q unitigs.hpc.fasta --pi 95 -s 10000 -f none -t 8 -o >(awk '$11 >= 50000' > unitigs.hpc.mashmap.out)
         $mashmap -r "$hpcRef" \
            -q unitigs.hpc.fasta \
            --pi 95 \
            -s 10000 \
            -f none \
            -t 8 \
            -o unitigs.hpc.mashmap.raw
      fi

      # awk '$11 >= 50000' unitigs.hpc.mashmap.raw > unitigs.hpc.mashmap.out
      awk '$11 >= 10000' unitigs.hpc.mashmap.raw > unitigs.hpc.mashmap.out
      cd ..
      ln -s 8-hicPipeline/unitigs.hpc.mashmap.out $utigs_mashmap
   else
      cd 5-untip
	  if [ ! -e unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.fasta ]; then
	     cat unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.gfa |awk '{if (match($1, "^S")) { print ">"$2; print $3}}'|fold -c  > unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.fasta
	  fi
	  # $mashmap -r $hpcRef -q unitig*.fasta --pi 95 -s 10000 -f none -t 8  -o >(awk '$11 >= 50000' > unitigs.hpc.mashmap.out)
     if [ ! -s unitigs.hpc.mashmap.raw ]; then
      $mashmap -r "$hpcRef" \
            -q unitig*.fasta \
            --pi 95 \
            -s 10000 \
            -f none \
            -t 8 \
            -o unitigs.hpc.mashmap.raw
      fi

      awk '$11 >= 50000' unitigs.hpc.mashmap.raw > unitigs.hpc.mashmap.out
      cd ..
      ln -s 5-untip/unitigs.hpc.mashmap.out $utigs_mashmap
   fi
fi

hpcGraph="assembly.homopolymer-compressed.noseq.gfa"
if [ ! -e $hpcGraph ]; then
   if [ $isHIC -gt 0 ]; then
      ln -s 8-hicPipeline/unitigs.hpc.noseq.gfa assembly.homopolymer-compressed.noseq.gfa
   else
      ln -s 5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.gfa assembly.homopolymer-compressed.noseq.gfa
   fi
fi

chrName=$4

hasNC=`cat assembly.mashmap.out |grep -c ${chrName}`
g="."
if [ $hasNC -gt 0 ]; then
   g=${chrName}
fi
hasCNC=`cat $utigs_mashmap |grep -c NC`
cg="."
if [ $hasCNC -gt 0 ]; then
   cg="NC"
fi
label1="sire_"
label2="dam_"

isMAT=`grep -c "sire_" $contigs`
if [ $isMAT -eq 0 ]; then
   isHF=`grep -c "h1tg" $contigs`
   if [ $isHF -eq 0 ]; then
      isUnpaired=`grep -c "contig-" $contigs`
      if [ $isUnpaired -eq 0 ]; then
         label1="haplotype1"
         label2="haplotype2"
      else
         label1="contig-"
         label2="none-ignore"
      fi
   else
      label1="h1tg"
      label2="h2tg"
   fi
fi
echo "$isMAT $label1 $label2 compNC: $cg regNC: $g"

num_contig=$(grep -c $label1 $contigs)
if [ $num_contig -eq 0 ]; then
   printf "Error: no contigs found with label $label1 in $contigs\n"
   exit 1
fi

minLen=$5
NUM=`cat assembly.homopolymer-compressed.chr.csv 2>/dev/null |wc -l`
if [ $NUM -le 1 ]; then
   echo -e "node\tchr" > assembly.homopolymer-compressed.chr.csv
   for i in `cat $utigs_mashmap |awk '{for (i = 1; i<=NF; i++) { if (match($i, "id:f:")) IDY=substr($i, 6, length($i)); } print $0"\t"IDY*100}' | awk '{if ($NF > IDY) print $6}'|sort |uniq`; do
      echo "Chr $i"
      cat $utigs_mashmap |awk '{for (i = 1; i<=NF; i++) { if (match($i, "id:f:")) IDY=substr($i, 6, length($i)); } print $0"\t"IDY*100}' |awk -v M=$minLen -v IDY=$idy '{if ($NF > IDY && $4-$3 > M) print $1"\t"$6"\t"$2}'|grep -w $i |sort -srnk3,3 |awk '{print $1"\t"$2}' |sort |uniq | grep "$cg" >> assembly.homopolymer-compressed.chr.csv
   done

   cat assembly.homopolymer-compressed.chr.csv |awk '{print $1}' |sort |uniq > tmp4
   $neigh assembly.homopolymer-compressed.noseq.gfa tmp.gfa -n tmp4 --drop-sequence -r 1000
   cat tmp.gfa |grep "^S" |awk '{print $2}' > tmp4
   #second pass to get missing chr w/shorter matches
   minLen=`echo $minLen |awk '{print $1/10}'`
   cat $utigs_mashmap |awk '{for (i = 1; i<=NF; i++) { if (match($i, "id:f:")) IDY=substr($i, 6, length($i)); } print $0"\t"IDY*100}' |grep -w -v -f tmp4 |awk -v M=$minLen -v IDY=$idy '{if ($NF > IDY && $4-$3 > M) print $1"\t"$6}' |sort |uniq | grep "$cg" >> assembly.homopolymer-compressed.chr.csv
   rm tmp.gfa tmp4
fi
minLen=1000000
cat assembly.mashmap.out |awk '{for (i = 1; i<=NF; i++) { if (match($i, "id:f:")) IDY=substr($i, 6, length($i)); } print $0"\t"IDY*100}' |grep "$label1" |grep $g |awk -v M=$minLen -v IDY=$idy '{if ($NF > IDY && $4-$3 > M) print $1"\t"$6"\t"$2"\t"$7}' |sort |uniq > translation_hap1
cat assembly.mashmap.out |awk '{for (i = 1; i<=NF; i++) { if (match($i, "id:f:")) IDY=substr($i, 6, length($i)); } print $0"\t"IDY*100}' |grep "$label2" |grep $g |awk -v M=$minLen -v IDY=$idy '{if ($NF > IDY && $4-$3 > M) print $1"\t"$6"\t"$2"\t"$7}' |sort |uniq > translation_hap2

minLen=`echo $minLen |awk '{print $1*15}'`
cat translation_hap1 | sort -k2,2|awk -v M=$minLen '{if ($3 > M) print $0}' |awk -v LAST="" -v S="" '{if (LAST != $2) { if (S > 0) print LAST"\t"C"\t"SUM/S*100"\t"MAX/S*100"\t"TIG; SUM=0; MAX=0; C=0; } LAST=$2; S=$NF; SUM+=$3; if (MAX < $3) MAX=$3; C+=1; TIG=$1} END {print LAST"\t"C"\t"SUM/S*100"\t"MAX/S*100"\t"TIG;}'  |awk '{print $1"\t"$4}' |sort -nk1,1 -s > chr_completeness_max_hap1

if [ -s translation_hap2 ]; then
   cat translation_hap2 | sort -k2,2|awk -v M=$minLen '{if ($3 > M) print $0}' |awk -v LAST="" -v S="" '{if (LAST != $2) { if (S > 0) print LAST"\t"C"\t"SUM/S*100"\t"MAX/S*100"\t"TIG; SUM=0; MAX=0; C=0; } LAST=$2; S=$NF; SUM+=$3; if (MAX < $3) MAX=$3; C+=1; TIG=$1} END {print LAST"\t"C"\t"SUM/S*100"\t"MAX/S*100"\t"TIG;}'  |awk '{print $1"\t"$4}' |sort -nk1,1 -s > chr_completeness_max_hap2
fi