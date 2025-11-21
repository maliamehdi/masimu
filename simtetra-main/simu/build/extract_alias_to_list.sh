#!/usr/bin/env bash
set -euo pipefail
# Extract a plain one-per-line energy list from an alias macro
# Usage: extract_alias_to_list.sh INPUT_ALIAS_MAC [OUTPUT_LIST]
# If OUTPUT_LIST is omitted, prints to stdout.

infile=${1:?Need input alias macro}
outfile=${2:-}

if [[ ! -s "$infile" ]]; then
  echo "Input not found or empty: $infile" >&2
  exit 1
fi

# Strategy: extract numbers in order (ints or floats)
if [[ -z "$outfile" ]]; then
  grep -Eo '[0-9]+(\.[0-9]+)?' "$infile"
else
  grep -Eo '[0-9]+(\.[0-9]+)?' "$infile" > "$outfile"
  echo "Wrote energy list to $outfile" >&2
fi
