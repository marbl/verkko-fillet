#!/bin/sh

module load bedtools/2.31.1

MIN_LEN=200000
if [[ "$#" -ge 2 ]]; then
  MIN_LEN=$2
fi
MIN_FRAC=10
if [[ "$#" -ge 3 ]]; then
   MIN_FRAC=$3
fi
MAPPED=""

IGNORED_BED="$tools/verkko-fillet/src/verkkofillet/data/dataset/chm13v2.0.ignore.bed"
MASHMAP=$1

# Extract unique query names that pass the length threshold
queries=$(awk -v M="$MIN_LEN" '$2 >= M {print $1}' "$MASHMAP" | sort -u)

for q in $queries; do
  #   1. Grab this query's mashmap rows, strip "id:f:" tags
  #   2. Keep high-identity hits (>0.99); prepend ref coords for bedtools
  #   3. Drop hits that fall entirely inside an ignored region
  prefiltered=$(
    grep -w "$q" "$MASHMAP" \
    | sed 's/id:f://g' \
    | awk '$(NF-1) > 0.99 {print $6"\t"$8"\t"$9"\t"$0}' \
    | bedtools intersect -f 1 -v -a - -b "$IGNORED_BED"
  )

  if [ -z "$prefiltered" ]; then
    echo "WARNING: no mapping >99% for $q" >&2
    continue
  fi

  #   4. Reshape to: query, qlen, qspan, strand, ref, rspan
  #   5. Sort by ref so we can aggregate per-ref in one pass
  #   6. For each ref, sum forward/reverse span and emit % coverage
  #   7. Decide orientation: if >MIN_FRAC covered, emit a BED-ish assignment line
  result=$(
    printf '%s\n' "$prefiltered" \
    | awk '{print $4"\t"$5"\t"$7-$6"\t"$8"\t"$9"\t"$12-$11}' \
    | sort -sk5,5 \
    | awk '
      function flush() {
        if (prev != "") print prev"\t"name"\t"fwd/len*100"\t"rev/len*100"\t"len
      }
      {
        if ($5 != prev) { flush(); prev=$5; name=$1; len=$2; fwd=0; rev=0 }
        if ($4 == "-") rev += $3; else fwd += $3
      }
      END { flush() }
    ' \
    | awk -v MF="$MIN_FRAC" '
      $3 > MF && $3 >  $4 { print $2"\t0\t"$NF"\t"$1"\t"$3"\t+" }
      $4 > MF && $4 >= $3 { print $2"\t0\t"$NF"\t"$1"\t"$4"\t-" }
    '
  )

  # Count non-empty assignment lines
  result=$(printf '%s\n' "$result" | grep -v '^[[:space:]]*$')
  n=$(printf '%s\n' "$result" | grep -c .)
  # echo $n $(echo "$result" | tr '\n' ' ') >&2

  if [ "$n" -eq 0 ]; then
    echo "WARNING: could not assign $q" >&2
  elif [ "$n" -gt 1 ]; then
    joined=$(printf '%s\n' "$result" | sort -sk5,5nr \
    | awk '{
      # First row sets the shared BED coords (col 1-3) and tracks max score
      if (NR == 1) { c1=$1; c2=$2; c3=$3; best=$5; dir=$6 }
      score += $5
      name = (name ? name"_" : "") $4
    }
    END { print c1"\t"c2"\t"c3"\t"name"_"c1"\t"score"\t"dir } ')
    echo "WARNING: more than one assignment for $q, it has $(echo "$result" | tr '\n' ' ')" >&2
    echo "$joined"
  else
    echo "$result" | awk '{print $1"\t"$2"\t"$3"\t"$4"_"$1"\t"$5"\t"$6}'
  fi
done
