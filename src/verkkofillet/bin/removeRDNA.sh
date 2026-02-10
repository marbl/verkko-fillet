#!/bin/bash
# will idenxify sequences matching user-provided repeat unit or rDNA by default
# will also mark and add telomere sequences on the graph
# requires $mash and seqtk along with verkko in path

mash=$(which mash 2>/dev/null)
if [ "x$mash" == "x" ]; then
   module load mash
   mash=$(which mash 2>/dev/null)
fi
if [ "x$mash" == "X" ]; then
   echo "Error: mash not found"
   exit
fi
seqtk=$(which seqtk 2>/dev/null)
if [ "x$seqtk" == "x" ]; then
   module load seqtk
   seqtk=$(which seqtk 2>/dev/null)
fi
if [ "x$seqtk" == "x" ]; then
   echo "Error: seqtk not found"
   exit
fi
verkko=$(which verkko 2>/dev/null)
if [ "x$verkko" == "x" ]; then
   verkko="/data/korens/devel/verkko-tip/bin/verkko"

   if [ ! -e $verkko ]; then
      echo "Error: verkko not found"
      exit
   fi
fi
verkko=`dirname $verkko |awk '{print $1"/../lib/verkko/scripts/remove_nodes_add_telomere.py"}'`
if [ ! -e $verkko ]; then
   echo "Error: no verkko src found"
   exit
fi

echo "Running with mash $mash seqtk $seqtk and verkko $verkko"
SCRIPT_DIR=$(dirname "$(realpath "${BASH_SOURCE[0]}")")

repeatUnit="$SCRIPT_DIR/../data/dataset/rDNA_compressed.fasta"
if [[ "$#" -ge 1 ]]; then
   echo "Using custom repeat sequence $1"
   $seqtk hpc $1 > repeatUnit.hpc.fasta
   repeatUnit="repeatUnit.hpc.fasta"
fi

isRUK=`ls 6-rukki/unitig*.fasta 2>/dev/null |wc -l |awk '{print $1}'`
isHIC=`ls 8-hicPipeline/unitigs.hpc.fasta 2> /dev/null |wc -l |awk '{print $1}'`
params=""

if [ ! -s compressed.sketch.msh ]; then
   if [ $isRUK -gt 0 ]; then
      cd 6-rukki
      $mash sketch -i unitig*.fasta -o sketch.msh
      cd ..
      ln -s 6-rukki/sketch.msh compressed.sketch.msh
   elif [ $isHIC -gt 0 ]; then
      cd 8-hicPipeline
      $mash sketch -i unitigs.hpc.fasta -o sketch.msh
      cd ..
      ln -s 8-hicPipeline/sketch.msh compressed.sketch.msh
   else
      cd 5-untip
      isDone=`ls unitig*.fasta 2>/dev/null |wc -l |awk '{print $1}'`
      if [ $isDone -le 0 ]; then
         cat unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.gfa |awk '{if (match($1, "^S")) { print ">"$2; print $3}}'|fold -c > unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.fasta
      fi
      $mash sketch -i unitig*.fasta -o sketch.msh
      cd ..
      ln -s 5-untip/sketch.msh compressed.sketch.msh
   fi
fi

if [ ! -e stats/assembly.telomere.bed ]; then
   echo "Need seqtk telo, please compute telomere first"
   exit
fi
$mash screen compressed.sketch.msh $repeatUnit | awk '{if ($1 > 0.9 && $4 < 0.05) print $NF}' > target.screennodes.out

python $verkko -r target.screennodes.out -t stats/assembly.telomere.bed
