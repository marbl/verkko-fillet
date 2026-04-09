#! /bin/bash

seqtk=$(which seqtk 2>/dev/null)
if [ "x$seqtk" == "x" ]; then
   module load seqtk
   seqtk=$(which seqtk 2>/dev/null)
fi
verkko=$(which verkko 2>/dev/null)
if [ "x$verkko" == "x" ]; then
   module load verkko
   verkko=$(which verkko 2>/dev/null)
fi


verkkofilletdir=$1
verkkodir=$2
newfolder=$3
finalGaf=$4
missing_edge_dir=$5

finalGaf=$(realpath $finalGaf)
newfolder=$(realpath $newfolder)
verkkofilletdir=$(realpath $verkkofilletdir)
verkkodir=$(realpath $verkkodir)
subsetid=$(realpath $missing_edge_dir)

echo -e "finalGaf : $finalGaf"
echo -e "newfolder : $newfolder"
echo -e "verkkofilletdir : $verkkofilletdir"
echo -e "verkkodir : $verkkodir"
echo -e "missing_edge_dir : $missing_edge_dir"


for dir in $newfolder $verkkodir $verkkofilletdir $missing_edge_dir; do
    if [ ! -d $dir ]; then
        echo "Error: $dir does not exist or is not a directory"
        exit 1
    fi
done

if [ ! -f $missing_edge_dir/ont_subset.tmp.fasta ]; then
    if [ -f $missing_edge_dir/*fa ]; then
        echo "fasta files are in the missing edge directory, concat them to create ont_subset.tmp.fasta"
        cat $missing_edge_dir/*fa > $missing_edge_dir/ont_subset.tmp.fasta
    else
        echo "No fasta files found in the missing edge directory, extracting fasta sequences from ONT data for missing edges..."
        echo "Extracting fasta sequences for missing edges..."
        cat $missing_edge_dir/patch.gapid*id > $missing_edge_dir/ont_subset.tmp.id
        cat $missing_edge_dir/patch.gapid*gaf > $missing_edge_dir/ont_subset.tmp.gaf
        zcat $verkkodir/3-align/split/ont*.fasta.gz | $seqtk subseq - $missing_edge_dir/ont_subset.tmp.id > $missing_edge_dir/ont_subset.tmp.fasta
    fi
else 
    echo "ont_subset.tmp.fasta already exists, skipping fasta extraction"
fi

echo "processing 7-consensus directory..."
mkdir -p $newfolder/7-consensus &&
cd $newfolder/7-consensus/ &&
cp $verkkodir/7-consensus/ont_subset.* ./ &&
chmod a+w * &&
cat $missing_edge_dir/*.missing_edge.ont_list.txt | uniq >> ont_subset.id &&
gunzip ont_subset.fasta.gz &&
cat $missing_edge_dir/ont_subset.tmp.fasta >> ont_subset.fasta &&
bgzip ont_subset.fasta &&
echo "7-consensus directory is updated"
cd ..

echo "processing 6-layoutContigs directory..."
cp -r $verkkodir/6-layoutContigs/ .
chmod -R a+w 6-layoutContigs/ &&
cd 6-layoutContigs/
chmod a+w * &&
rm consensus_paths.txt &&
cat $verkkofilletdir/missing_edge/patch.gapid_*.gaf >> hifi.alignments.gaf &&
cat $verkkofilletdir/missing_edge/gapid_*.missing_edge.patching.gfa | grep "^L" | grep gap >> combined-edges.gfa &&
cat $verkkofilletdir/missing_edge/gapid_*.missing_edge.patching.gfa | grep gap | awk 'BEGIN { FS="[ \t]+"; OFS="\t"; } ($1 == "S") && ($3 != "*") { print $2, length($3); }' >> nodelens.txt &&
cp $finalGaf ./consensus_paths.txt &&
cat $missing_edge_dir/ont_subset.tmp.id >> ont-gapfill.txt &&
rm ./unitig*

echo " "
echo "running replace_path_nodes.py"
echo " "
verkkoLib=$($verkko | grep "Verkko module"| awk '{print $4}')
$verkkoLib/scripts/replace_path_nodes.py ../4-processONT/alns-ont-mapqfilter.gaf combined-nodemap.txt |grep -F -v -w -f ont-gapfill.txt > ont.alignments.gaf &&
cd ..

echo "6-layoutContigs directory is updated!"
