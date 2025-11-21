#!/usr/bin/env bash
set -euo pipefail
# Usage: gen_macro_range.sh MIN_KEV MAX_KEV STEP_KEV BEAMON [OUTFILE]
# Example: ./gen_macro_range.sh 30 140000 30 5000 prompt_PARIS50_30to140000_30keV.mac

if [[ ${1:-} == "-h" || $# -lt 4 ]]; then
  echo "Usage: $0 MIN_KEV MAX_KEV STEP_KEV BEAMON [OUTFILE]" >&2
  exit 1
fi

MIN=$1
MAX=$2
STEP=$3
BEAM=$4
OUT=${5:-/dev/stdout}

# Write to temp then move atomically if OUT is a file
if [[ "$OUT" != "/dev/stdout" ]]; then
  TMP="$OUT.tmp.$$"
else
  TMP="/dev/stdout"
fi

{
  echo "# Auto-generated Geant4 macro: ${MIN}keV to ${MAX}keV step ${STEP}keV"
  echo "/control/verbose 0"
  echo "/run/verbose 0"
  echo "/event/verbose 0"
  echo "/tracking/verbose 0"
  echo "/run/initialize"
  echo "/gun/particle gamma"
  echo "/gun/position 0 0 0 mm"
  echo "/gun/direction 0 0 1"
  
  # Main loop (includes last multiple <= MAX)
  for e in $(seq "$MIN" "$STEP" "$MAX"); do
    echo "/gun/energy ${e} keV"
    echo "/run/beamOn ${BEAM}"
    echo
  done
  
  # If MAX is not a multiple of STEP, optionally append exact MAX as last point
  if (( MAX % STEP != 0 )); then
    echo "# Append exact MAX not on STEP grid"
    echo "/gun/energy ${MAX} keV"
    echo "/run/beamOn ${BEAM}"
  fi
} > "$TMP"

if [[ "$OUT" != "/dev/stdout" ]]; then
  mv "$TMP" "$OUT"
  echo "Wrote macro to $OUT" >&2
fi
