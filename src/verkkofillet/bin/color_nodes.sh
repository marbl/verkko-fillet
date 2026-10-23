#!/bin/sh

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 node_ids.txt label color"
  echo "  node_ids.txt  Path to a file containing node IDs. One node ID per line at the last column."
  echo "  label         Label to assign to the nodes."
  echo "  color         Color to assign to the nodes. e.g. \'#A020F0\' for rDNA, \'#008000\' for telomere."
  exit 1
fi

node_ids_file=$1
label=$2
color=$3

if [ ! -f "$node_ids_file" ]; then
  echo "Error: Node IDs file $node_ids_file does not exist."
  exit 1
fi

echo -e "node\t$label\tcolor" > $node_ids_file.$label.csv
awk -v label=$label -v color=$color '{print $NF"\t"label"\t"color}' $node_ids_file >> $node_ids_file.$label.csv
