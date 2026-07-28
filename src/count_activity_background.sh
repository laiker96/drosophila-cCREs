#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 4 ]; then
  echo "usage: count_activity_background.sh UNITS.bed.gz BINS.bed OUTPUT.tsv.gz THREADS" >&2
  exit 2
fi

units=$1
bins=$2
output=$3
threads=$4
mkdir -p "$(dirname "$output")"
temporary=$(mktemp -d "$(dirname "$output")/.$(basename "$output").XXXXXX")
trap 'rm -rf "$temporary"' EXIT

pigz -dc "$units" \
  | awk 'BEGIN {OFS="\t"}
      NR==FNR {
        n++; chrom[n]=$1; start[n]=$2; end[n]=$3; id[n]=$4;
        bin_index[$1 SUBSEP int($2/10000)]=n; next
      }
      {
        middle=int(($2+$3)/2); key=$1 SUBSEP int(middle/10000);
        if (key in bin_index) count[bin_index[key]]++
      }
      END {
        print "background_bin_id","chrom","start","end","raw_count";
        for (i=1; i<=n; i++) print id[i],chrom[i],start[i],end[i],count[i]+0
      }' "$bins" - \
  | pigz -p "$threads" -c > "$temporary/counts.tsv.gz"

pigz -t "$temporary/counts.tsv.gz"
test "$(pigz -dc "$temporary/counts.tsv.gz" | wc -l)" \
  -eq "$(( $(wc -l < "$bins") + 1 ))"
mv "$temporary/counts.tsv.gz" "$output"
