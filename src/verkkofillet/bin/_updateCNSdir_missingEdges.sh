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
if [ "x$samtools" == "x" ]; then
   module load samtools
   samtools=$(which samtools 2>/dev/null)
fi


verkkofilletdir=$1
verkkodir=$2
newfolder=$3
finalGaf=$4
missing_edge_dir=$5
verkko_version=$6

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
echo -e "verkko_version : $verkko_version"

## Clean up any existing temporary files in the missing edge directory
rm -rf $missing_edge_dir/ont_subset.tmp.id $missing_edge_dir/ont_subset.tmp.gaf $missing_edge_dir/ont_subset.tmp.fasta
error_count=0

## Check if the required directories and files exist
for dir in $newfolder $verkkodir $verkkofilletdir $missing_edge_dir; do
    if [ ! -d $dir ]; then
        echo "Error: $dir does not exist or is not a directory"
        exit 1
    fi
done

## Generate the list of ONT read IDs and corresponding GAF entries for the missing edges
cat $missing_edge_dir/*.missing_edge.ont_list.txt | sort | uniq > $missing_edge_dir/ont_subset.tmp.id
cat $missing_edge_dir/patch.gapid*gaf > $missing_edge_dir/ont_subset.tmp.gaf

## Generate the fasta file for the ONT reads in the missing edges
## Generate the fasta file for the ONT reads in the missing edges
rm -rf $missing_edge_dir/ont_subset.tmp.fasta $missing_edge_dir/ont_subset.tmp.fasta.fai

fa_count=$(find "$missing_edge_dir" -maxdepth 1 -name '*.fa' 2>/dev/null | wc -l)
fasta_count=$(find "$missing_edge_dir" -maxdepth 1 -name '*.fasta' ! -name 'ont_subset.tmp.fasta' 2>/dev/null | wc -l)

if [ "$fa_count" -gt 0 ] || [ "$fasta_count" -gt 0 ]; then
    echo "fasta files are in the missing edge directory, concat them to create ont_subset.tmp.fasta"
        find "$missing_edge_dir" -maxdepth 1 \( -name '*.fa' -o -name '*.fasta' \) ! -name 'ont_subset.tmp.fasta'
    find "$missing_edge_dir" -maxdepth 1 \( -name '*.fa' -o -name '*.fasta' \) ! -name 'ont_subset.tmp.fasta' -exec cat {} + > "$missing_edge_dir/ont_subset.tmp.fasta"
    samtools faidx $missing_edge_dir/ont_subset.tmp.fasta
fi
cut -f 1 $missing_edge_dir/ont_subset.tmp.fasta.fai > $missing_edge_dir/ont_subset.tmp.fasta.fai.id
grep -w -v -f $missing_edge_dir/ont_subset.tmp.fasta.fai.id $missing_edge_dir/ont_subset.tmp.id > $missing_edge_dir/ont_subset.tmp.id.missing
rest_read_count=$(wc -l < $missing_edge_dir/ont_subset.tmp.id.missing)

if [ "$rest_read_count" -gt 0 ]; then
    echo "$rest_read_count ONT reads in ont_subset.tmp.id are not found in the fasta files in the missing edge directory, extracting them from the original fasta files in 3-align/split/ont*.fasta.gz"
    # zcat "$verkkodir"/3-align/split/ont*.fasta.gz | "$seqtk" subseq - "$missing_edge_dir/ont_subset.tmp.id.missing" >> "$missing_edge_dir/ont_subset.tmp.fasta"
    find "$verkkodir"/3-align/split/ -name 'ont*.fasta.gz' | \
    parallel -j 8 "zcat {} | $seqtk subseq - $missing_edge_dir/ont_subset.tmp.id.missing" \
    >> "$missing_edge_dir/ont_subset.tmp.fasta"
fi
rm -rf $missing_edge_dir/ont_subset.tmp.id.missing $missing_edge_dir/ont_subset.tmp.fasta.fai $missing_edge_dir/ont_subset.tmp.fasta.fai.id
samtools faidx $missing_edge_dir/ont_subset.tmp.fasta

## Check if the generated files are not empty
check_files="$missing_edge_dir/ont_subset.tmp.id $missing_edge_dir/ont_subset.tmp.gaf $missing_edge_dir/ont_subset.tmp.fasta"
for file in $check_files; do
    if [ ! -s $file ]; then
        echo "Error: $file is empty, please check the corresponding step for generating this file"
        exit 1
    fi
done

## Run for 7-consensus
echo "processing 7-consensus directory..."
mkdir -p $newfolder/7-consensus &&
cd $newfolder/7-consensus/ &&
cp $verkkodir/7-consensus/ont_subset.* ./ &&
chmod a+w * &&
cat $missing_edge_dir/ont_subset.tmp.id | uniq >> ont_subset.id &&
gunzip ont_subset.fasta.gz &&
cat $missing_edge_dir/ont_subset.tmp.fasta >> ont_subset.fasta &&
bgzip ont_subset.fasta &&
samtools faidx ont_subset.fasta.gz &&
echo "7-consensus directory is updated"
cd ..


## Run for 6-layoutContigs
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

## Run replace_path_nodes.py to update the consensus_paths.txt with the new nodes from the missing edges, and filter out the ont alignments that are not in ont-gapfill.txt
echo " "
echo "running replace_path_nodes.py"
echo " "
verkkoLib=$($verkko | grep "Verkko module"| awk '{print $4}')
$verkkoLib/scripts/replace_path_nodes.py ../4-processONT/alns-ont-mapqfilter.gaf combined-nodemap.txt |grep -F -v -w -f ont-gapfill.txt > ont.alignments.gaf  || exit 1
cd ..

# Check if all ONT reads in ont-gapfill.txt are included in output files

# check read id
files_to_check=(
    6-layoutContigs/hifi.alignments.gaf
    7-consensus/ont_subset.id
    7-consensus/ont_subset.fasta.gz.fai
    6-layoutContigs/ont-gapfill.txt
)

for file in "${files_to_check[@]}"; do
    if [ ! -s $file ]; then
        echo "Error: $file is empty, please check the corresponding step for generating this file"
        # exit 1
        error_count=$((error_count + 1))
    else
        # check if all ONT reads (ont_list) in file
        for ont_read in $(cat $missing_edge_dir/ont_subset.tmp.id); do
            if ! grep -qw "$ont_read" "$file"; then
                echo "Error: ONT read ID $ont_read is missing in $file"
                # exit 1
                error_count=$((error_count + 1))
            else
                echo -e "Checked: ONT read ID $ont_read is present in $file"
            fi
        done
    fi
done

# node id
files_to_check=(
    "6-layoutContigs/combined-edges.gfa"
    "6-layoutContigs/nodelens.txt"
)

cat $missing_edge_dir/gapid_*patching.gfa | cut -f 2 | grep "^gap" | sort | uniq  > "$missing_edge_dir/ont_subset.tmp.nodeid"
node_file="$missing_edge_dir/ont_subset.tmp.nodeid"
line_num_id=$(wc -l < "$node_file")

for file in "${files_to_check[@]}"; do
    if [ ! -s "$file" ]; then
        echo "Error: $file is empty. Please check the step that generates this file."
        # exit 1
        error_count=$((error_count + 1))
    else
        while read -r new_node; do
            if ! grep -qw "$new_node" "$file"; then
                echo "Error: node ID $new_node is missing in $file"
                # exit 1
                error_count=$((error_count + 1))
            else
                echo -e "Checked: $file contains all node IDs from $node_file."
            fi
        done < "$node_file"
        
    fi
done

for new_new in $(cat $missing_edge_dir/ont_subset.tmp.nodeid)
do
    if ! grep -q -w $new_new 6-layoutContigs/consensus_paths.txt; then
        echo "Error: node id $new_new is missing in 6-layoutContigs/consensus_paths.txt"
        # exit 1
        error_count=$((error_count + 1))
    fi
    echo -e "Checked: node id $new_new is present in 6-layoutContigs/consensus_paths.txt"
done

## check if the folders are not empty
folders_to_check=(
    1-buildGraph
    2-processGraph
    3-align
    4-processONT
    5-untip
    6-layoutContigs
    7-consensus
)

for folder in "${folders_to_check[@]}"; do
    if [ -d $folder ]; then
        if [ "$(ls -A $folder)" ]; then
            echo "Checked: $folder is not empty"
        else
            echo "Error: $folder is empty, please check the corresponding step for generating files in this folder"
            #exit 1
            error_count=$((error_count + 1))
        fi
    else
        echo "Error: $folder does not exist, please check the corresponding step for creating this folder"
        #exit 1
        error_count=$((error_count + 1))
    fi
done

echo -e "Total errors found: $error_count"
if [ $error_count -gt 0 ]; then
    echo "Please address the above errors before proceeding to the next steps."
else
    echo "All checks passed successfully. You can proceed to the next steps."
fi