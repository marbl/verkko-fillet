#!/bin/sh

module load mashmap/3.1.1
module load mash/2.3
module load samtools/1.23
module load python/3.10
module load bedtools/2.31.1
module load seqtk/1.5

cpus=$SLURM_CPUS_PER_TASK
if [[ -z $cpus ]]; then
  cpus=24
fi

tools="/data/Phillippy/tools"
refn="$tools/verkko-fillet/src/verkkofillet/data/dataset/chm13v2.0.fa.gz"
refc="$tools/verkko-fillet/src/verkkofillet/data/dataset/chm13v2.0.hpc.fa.gz"
rdna="$tools/verkko-fillet/src/verkkofillet/data/dataset/rDNA_compressed.fasta"

if [[ ! -s $refn ]]; then
  echo "Downloading $refn"
  wget -O $refn https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13v2.0.fa.gz
  wget -O $refn.gzi https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13v2.0.fa.gz.gzi
  wget -O $refn.fai https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13v2.0.fa.gz.fai
fi
if [[ ! -s $refc ]]; then
  echo "Downloading $refc"
  wget -O $refc https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13v2.0.hpc.fa.gz
  wget -O $refc.gzi https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13v2.0.hpc.fa.gz.gzi
  wget -O $refc.fai https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13v2.0.hpc.fa.gz.fai
fi

mkdir -p analysis && cd analysis

if [ ! -f ../assembly.fasta.fai ]; then
  echo "Index ../assembly.fasta"
  samtools faidx -@${cpus} ../assembly.fasta
fi

# collect telomeres, gaps and T2T status
if [ ! -e assembly.t2t_ctgs ]; then
  sh $tools/verkko-fillet/src/verkkofillet/bin/getT2T.sh ../assembly.fasta
fi

# remove rDNA and add telomere nodes
if [ ! -s assembly.homopolymer-compressed.add_telo.noseq.gfa ]; then
   ln -s ../assembly.colors.csv
   cat ../8-hicPipeline/unitigs.hpc.fasta | mash sketch -i - -o sketch.msh
   mash screen sketch.msh $rdna | awk '{if ($1 > 0.9 && $4 < 0.05) print $NF}' > rdna.nodes
   python $tools/verkko-fillet/src/verkkofillet/bin/remove_nodes_add_telomere.py \
     -t assembly.telomere.bed \
     -g ../assembly.homopolymer-compressed.noseq.gfa \
     -s ../assembly.scfmap -p ../assembly.paths.tsv \
     -o assembly.homopolymer-compressed.add_telo.noseq.gfa  \
     -c assembly.colors.add_telo_add_rdna.csv
   python $tools/verkko-fillet/src/verkkofillet/bin/remove_nodes_add_telomere.py \
     -t assembly.telomere.bed \
     -r rdna.nodes \
     -g ../assembly.homopolymer-compressed.noseq.gfa \
     -s ../assembly.scfmap -p ../assembly.paths.tsv \
     -o assembly.homopolymer-compressed.add_telo_remove_rdna.noseq.gfa \
     -c assembly.colors.add_telo_add_rdna.csv
   unlink ./assembly.colors.csv
fi

if [ ! -e assembly-ref.norm.mashmap ]; then
  mashmap -r $refn -q ../assembly.fasta --pi 95 -s 10000 -t ${cpus} -o assembly-ref.norm.mashmap
fi

if [ ! -e assembly-ref.comp.mashmap ]; then
  mashmap -r $refc -q ../8-hicPipeline/unitigs.hpc.fasta --pi 95 -s 10000 -t ${cpus} -o assembly-ref.comp.mashmap
fi

g="."
if grep -q NC assembly-ref.norm.mashmap ; then
  g="NC"
fi

cg="."
if grep -q NC assembly-ref.comp.mashmap ; then
  cg="CNC"
fi

if grep -q mat- ../assembly.fasta ; then
  label1="pat-"
  label2="mat-"
elif grep -q h1tg ../assembly.fasta ; then
  label1="h1tg"
  label2="h2tg"
elif grep -q contig- ../assembly.fasta ; then
  label1="contig-"
  label2="none-ignore"
elif grep -q "#1#" ../assembly.fasta ; then
   label1="#1#"
   label2="#2#"
else
  label1="haplotype1-"
  label2="haplotype2-"
fi
echo "$label1 $label2 compNC: $cg regNC: $g"

# If nc-to-chr.sed is available, use that instead. Otherwise an empty nc-to-chr.sed will be used
if [[ ! -s nc-to-chr.sed ]]; then
  echo  > nc-to-chr.sed
fi

if [ ! -e assembly.homopolymer-compressed.chr.csv ] ; then
  echo -e "node\tchr" > assembly.homopolymer-compressed.chr.csv
  for i in $( cat assembly-ref.comp.mashmap | sed 's/id:f://g' | awk '{idy=$(NF-1)*100; if (idy > 99) print $6}'| sort -Vu ); do
    echo "Chr $i"
    cat assembly-ref.comp.mashmap | sed 's/id:f://g' | \
      awk '{idy=$(NF-1)*100; if (idy > 99 && $4-$3 > 5000000) print $1"\t"$6"\t"$2}' | \
      grep -w $i |sort -srnk3,3 | \
      awk '{print $1"\t"$2}' | \
      sort -u | \
      grep "$cg" | \
      sed -f nc-to-chr.sed >> assembly.homopolymer-compressed.chr.csv
  done

  cat assembly.homopolymer-compressed.chr.csv | awk '{print $1}' | sort -u > tmp

  python $tools/verkko-fillet/src/verkkofillet/bin/neighborhood.py -n tmp -r 1000 ../assembly.homopolymer-compressed.noseq.gfa > tmp.neighborhood

  #second pass to get missing chr w/shorter matches
  cat assembly-ref.comp.mashmap | sed 's/id:f://g' |\
    grep -w -v -f tmp.neighborhood | \
    awk '{idy=$(NF-1)*100; if (idy > 99 && $4-$3 > 500000) print $1"\t"$6}' | \
    sort -Vu | \
    grep "$cg" | \
    sed -f nc-to-chr.sed >> assembly.homopolymer-compressed.chr.csv

  # clean up
  rm tmp tmp.neighborhood
fi

if [ ! -e translation_hap1 -o ! -e translation_hap2 ] ; then
  cat assembly-ref.norm.mashmap | grep "$label1" | grep $g \
    | sed 's/id:f://g' | awk '{idy=$(NF-1)*100; if (idy > 99 && $4-$3 > 1000000) print $1"\t"$6"\t"$2"\t"$7}' \
    | sort -Vu | sed -f nc-to-chr.sed > translation_hap1
  cat assembly-ref.norm.mashmap | grep "$label2" | grep $g \
    | sed 's/id:f://g' | awk '{idy=$(NF-1)*100; if (idy > 99 && $4-$3 > 1000000) print $1"\t"$6"\t"$2"\t"$7}' \
    | sort -Vu | sed -f nc-to-chr.sed > translation_hap2
fi

if [ ! -e chr_completeness_max_hap1 -o ! -e chr_completeness_max_hap2 ] ; then
  cat translation_hap1 |sort -k2,2V | awk '{if ($3 > 15000000) print $0}' \
    | awk -v LAST="" -v S="" '{if (LAST != $2) 
                                 { if (S > 0) print LAST"\t"C"\t"SUM/S*100"\t"MAX/S*100"\t"TIG;
                                   SUM=0; MAX=0; C=0; }
                              LAST=$2; S=$NF; SUM+=$3; if (MAX < $3) MAX=$3; C+=1; TIG=$1
                              } END {print LAST"\t"C"\t"SUM/S*100"\t"MAX/S*100"\t"TIG;}' \
    | awk '{print $1"\t"$4}' | sort -nk1,1 -s > chr_completeness_max_hap1
  cat translation_hap2 |sort -k2,2V | awk '{if ($3 > 15000000) print $0}' \
    | awk -v LAST="" -v S="" '{if (LAST != $2) 
                                 { if (S > 0) print LAST"\t"C"\t"SUM/S*100"\t"MAX/S*100"\t"TIG;
                                   SUM=0; MAX=0; C=0; }
                              LAST=$2; S=$NF; SUM+=$3; if (MAX < $3) MAX=$3; C+=1; TIG=$1
                              } END {print LAST"\t"C"\t"SUM/S*100"\t"MAX/S*100"\t"TIG;}' \
    | awk '{print $1"\t"$4}' | sort -nk1,1 -s > chr_completeness_max_hap2
fi

if [ ! -e assembly.refOriented.fasta.gz ]; then
  set -x
  sh $tools/verkko-fillet/src/verkkofillet/bin/reorientByRef.sh assembly-ref.norm.mashmap > assembly-ref.reorient.bed
  # bedtools getfasta -s adds (+) or (-) to the sequence headers. let's avoid this
  cat assembly-ref.reorient.bed | awk '$NF=="+" {print $1}' > assembly-ref.reorient.fwd
  cat assembly-ref.reorient.bed | awk '$NF=="-" {print $1}' > assembly-ref.reorient.rev
  samtools faidx -@$cpus -r    assembly-ref.reorient.fwd ../assembly.fasta >  tmp.fasta
  samtools faidx -@$cpus -i --mark-strand no -r assembly-ref.reorient.rev ../assembly.fasta >> tmp.fasta
  bedtools getfasta -fi tmp.fasta -bed assembly-ref.reorient.bed -nameOnly -fo assembly.refOriented.fasta
  bgzip -@$cpus assembly.refOriented.fasta
  samtools faidx -@$cpus assembly.refOriented.fasta.gz
  grep $label1 assembly.refOriented.fasta.gz.fai | cut -f1 > ${label1/-/}.chr.list
  grep $label2 assembly.refOriented.fasta.gz.fai | cut -f1 > ${label2/-/}.chr.list
  samtools faidx -@$cpus -r ${label1/-/}.chr.list -o assembly.refOriented.haplotype1.fasta.gz --write-index assembly.refOriented.fasta.gz
  samtools faidx -@$cpus -r ${label2/-/}.chr.list -o assembly.refOriented.haplotype2.fasta.gz --write-index assembly.refOriented.fasta.gz
  rm assembly-ref.reorient.fwd assembly-ref.reorient.rev tmp.fasta ${label1/-/}.chr.list ${label2/-/}.chr.list
fi

