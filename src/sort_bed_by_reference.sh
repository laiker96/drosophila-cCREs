#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: sort_bed_by_reference.sh CHROM_SIZES TEMP_DIR" >&2
    exit 2
fi

chrom_sizes=$1
temporary=$2

if [[ ! -r "$chrom_sizes" ]]; then
    echo "chromosome sizes file is not readable: $chrom_sizes" >&2
    exit 2
fi

mkdir -p "$temporary"

# Prefix each record with its reference-order chromosome rank. GNU sort then
# spills runs to TEMP_DIR while keeping its in-memory buffer bounded. The rank
# is removed before emitting the original BED record.
awk 'BEGIN { OFS="\t" }
     FILENAME == ARGV[1] {
         if (NF != 2 || $2 !~ /^[0-9]+$/ || $2 <= 0 || $1 in rank) {
             print "invalid chromosome sizes record at line " FNR > "/dev/stderr"
             failed=1
             exit 1
         }
         rank[$1]=++chromosome_count
         next
     }
     NF < 3 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $2 > $3 {
         print "invalid BED record at input line " FNR > "/dev/stderr"
         failed=1
         exit 1
     }
     !($1 in rank) {
         print "BED record uses unknown chromosome: " $1 > "/dev/stderr"
         failed=1
         exit 1
     }
     { print rank[$1],$0 }
     END {
         if (!failed && chromosome_count == 0) {
             print "chromosome sizes file is empty" > "/dev/stderr"
             exit 1
         }
     }' "$chrom_sizes" - \
    | LC_ALL=C sort --temporary-directory="$temporary" --buffer-size=512M \
        --parallel=2 --stable -k1,1n -k3,3n -k4,4n \
    | cut -f2-
