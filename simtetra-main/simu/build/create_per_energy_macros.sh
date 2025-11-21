#!/usr/bin/env bash
set -euo pipefail
# Create one Geant4 macro per energy from an alias-list .mac
# Usage:
#   ./create_per_energy_macros.sh PARIS50 ../bincenters_out_res/energies_list_PARIS50.mac 5000 [OUTDIR]
# Notes:
# - Extracts all numeric values from the file (works with '/control/alias Elist { ... }')
# - Macro filenames: <PARIS_ID>_E<energy>keV.mac ('.' replaced by 'p' for safety)
# - Each macro contains a single /run/beamOn block for that energy

if [[ ${1:-} == "-h" || $# -lt 3 ]]; then
  echo "Usage: $0 PARIS_ID ALIAS_FILE EVENTS [OUTDIR]" >&2
  exit 1
fi

PARIS_ID="$1"
ALIAS_FILE="$2"
EVENTS="$3"
OUTDIR="${4:-mac_per_energy}/${PARIS_ID}"

mkdir -p "$OUTDIR"

if [[ ! -s "$ALIAS_FILE" ]]; then
  echo "Input not found or empty: $ALIAS_FILE" >&2
  exit 2
fi

# Extract energies: pull all numbers (int or float) in order
mapfile -t ENERGIES < <(grep -Eo '[0-9]+(\.[0-9]+)?' "$ALIAS_FILE")

if (( ${#ENERGIES[@]} == 0 )); then
  echo "No energies found in $ALIAS_FILE" >&2
  exit 3
fi

echo "Generating ${#ENERGIES[@]} macros into $OUTDIR ..." >&2

count=0
for E in "${ENERGIES[@]}"; do
  # safe filename: replace '.' by 'p'
  ESAFE=${E//./p}
  MAC="$OUTDIR/${PARIS_ID}_E${ESAFE}keV.mac"
  cat > "$MAC" <<EOF
/control/verbose 0
/run/verbose 0
/event/verbose 0
/tracking/verbose 0

/run/initialize
/gun/particle gamma
/gun/position 0 0 -31.8 mm
/gun/direction 0 0 1
/gun/energy ${E} keV
/run/beamOn ${EVENTS}
EOF
  ((count++))
  # light progress
  if (( count % 200 == 0 )); then echo "  ... ${count} written" >&2; fi
done

echo "Done: ${count} macro files." >&2
