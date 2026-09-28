#!/bin/bash
# given an assembly, this script will produce a count of t2t contigs and scaffolds

asm=assembly.fasta
if [[ "$#" -ge 1 ]]; then
   asm=$1
fi
telo="CCCTAA"
if [[ "$#" -ge 2 ]]; then
   telo=$2
fi
PREFIX=`basename $asm |sed s/.fasta.gz//g          | sed s/.fasta//g          | sed s/.fa/.gaps.bed/g`
gaps=`basename $asm   |sed s/.fasta.gz/.gaps.bed/g | sed s/.fasta/.gaps.bed/g | sed s/.fa/.gaps.bed/g`
fai="$asm.fai"

echo "Processing $asm with gaps $gaps and fai $fai"
seqtk=$(which seqtk 2>/dev/null)
if [ "x$seqtk" == "x" ]; then
   module load seqtk
   seqtk=$(which seqtk 2>/dev/null)
fi
if [ "x$seqtk" == "x" ]; then
   echo "Error: seqtk not found"
   exit 1
fi

samtools=$(which samtools 2>/dev/null)
if [ "x$samtools" == "x" ]; then
   module load samtools
   samtools=$(which samtools 2>/dev/null)
fi
if [ "x$samtools" == "x" ]; then
   echo "Error: samtools not found"
   exit 1
fi

bedtools=$(which bedtools 2>/dev/null)
if [ "x$bedtools" == "x" ]; then
   module load bedtools
   bedtools=$(which bedtools 2>/dev/null)
fi
if [ "x$bedtools" == "x" ]; then
   echo "Error: bedtools not found"
   exit 1
fi

if [ ! -e $PREFIX.telomere.bed ] || [ $PREFIX.telomere.bed -ot $asm ] ; then
echo "Need telomere"
   $seqtk telo -m $telo -d 50000 $asm > $PREFIX.tmp
   # also trim off the start/end of the sequence and find telomere again because seqtk sometimes misses telomere in these cases
   $seqtk trimfq -b 2000 -e 2000 $asm > $PREFIX.tmp.fa
   $seqtk telo -m $telo -d 50000 $PREFIX.tmp.fa >> $PREFIX.tmp
   $seqtk trimfq -b 4000 -e 4000 $asm > $PREFIX.tmp.fa
   $seqtk telo -m $telo -d 50000 $PREFIX.tmp.fa >> $PREFIX.tmp
   cat $PREFIX.tmp | $bedtools sort |$bedtools merge -c 4 -o max > $PREFIX.telomere.bed
fi

if [ ! -e $gaps ] || [ $gaps -ot $asm ]; then
echo "Need gaps"
  $seqtk gap $asm > $gaps
fi

if [ ! -e $fai ] || [ $fai -ot $asm ]; then
   $samtools faidx $asm
fi

# get list of nodes with gaps
cat $gaps|awk '{print $1}'|sort |uniq > $PREFIX.tmp
# exclude those nodes from the full list of scaffolds in the asm
grep -w -v -f $PREFIX.tmp $fai |awk '{print $1}' > $PREFIX.tmp2

# get list of things with >1 telomere
cat $PREFIX.telomere.bed |awk '{print $1}'|sort |uniq -c|awk '{if ($1 > 1) print $NF}' > $PREFIX.tmp

# get list of nodes w/>1 telomere and no gaps
grep -w -f $PREFIX.tmp $PREFIX.tmp2 > $PREFIX.tmp3

# keep only those that have two telomeres near the ends
grep -w -f $PREFIX.tmp $PREFIX.telomere.bed |awk '{print $1}' > $PREFIX.tmp4
grep -w -f $PREFIX.tmp4 $asm.fai |awk '{print $1"\t"$2}' |sort -sk1,1 > $PREFIX.tmp5
grep -w -f $PREFIX.tmp $PREFIX.telomere.bed |sort -sk1,1 > $PREFIX.tmp4
join $PREFIX.tmp4 $PREFIX.tmp5 |awk -v PREV="" '{if ($1 != PREV) { if (PREV != "" && C >=1 && E >=1 ) print PREV; PREV=$1; C=0; E=0 } if ($2 < 10000) C++; if ($3 + 10000 > $NF && $NF > 50000) { E++; } } END { if (C >=1 && E >= 1) print PREV}' > $PREFIX.t2t_scfs

grep -w -f $PREFIX.tmp3 $PREFIX.telomere.bed |awk '{print $1}' > $PREFIX.tmp4
grep -w -f $PREFIX.tmp4 $asm.fai |awk '{print $1"\t"$2}' |sort -sk1,1 > $PREFIX.tmp5
grep -w -f $PREFIX.tmp3 $PREFIX.telomere.bed | sort -sk1,1 > $PREFIX.tmp4
join $PREFIX.tmp4 $PREFIX.tmp5 |awk -v PREV="" '{if ($1 != PREV) { if (PREV != "" && C >=1 && E >=1 ) print PREV; PREV=$1; C=0; E=0 } if ($2 < 10000) C++; if ($3 + 10000 > $NF && $NF > 50000 ) { E++; } } END { if (C >=1 && E >= 1) print PREV}' > $PREFIX.t2t_ctgs

N=`wc -l $PREFIX.t2t_ctgs |awk '{print $1}'`
if [ -s $PREFIX.t2t_ctgs ] ; then
   echo "Ungapped two telomere: $N"
   grep -w -f $PREFIX.t2t_ctgs $PREFIX.telomere.bed
else
    echo Ungapped two telomere: NONE
fi

N=`wc -l $PREFIX.t2t_scfs |awk '{print $1}'`
echo ""
if [ -s $PREFIX.t2t_scfs ]; then
   echo "Gapped two telomere: $N"
   grep -w -f $PREFIX.t2t_scfs $PREFIX.telomere.bed
else
    echo Gapped two telomere: NONE
fi
rm $PREFIX.tmp*
