#!/usr/bin/env bash
set -euo pipefail
# Run all available PARIS alias macros using the batch runner (few files).
# It extracts energies, runs one energy per process (temporary macro deleted),
# and stores outputs under ../../myanalyse/<PARIS_ID>/output_<PARIS_ID>_E<energy>keV.root
#
# Usage:
#   ./run_all_paris_energy_lists.sh EVENTS MAX_PROCS [ALIAS_DIR] [PARIS_IDS...]
# Examples:
#   ./run_all_paris_energy_lists.sh 5000 8                     # auto-scan bincenters_out_res
#   ./run_all_paris_energy_lists.sh 2000 16 bincenters_out_res PARIS50 PARIS90
#
# Env:
#   QUIET=1         suppress app logs
#   G4APP=./simTetra
#
# Notes:
# - Per-PARIS runs are executed sequentially; intra-PARIS parallelism is controlled by MAX_PROCS.
#   This avoids overcommitting cores by nesting parallelism.

EVENTS=${1:?Need EVENTS per energy}
MAX_PROCS=${2:?Need MAX_PROCS}
ALIAS_DIR=${3:-../bincenters_out_res}
shift 0 || true

here=$(cd "$(dirname "$0")" && pwd)
runner="$here/run_energy_list.sh"

if [[ ! -x "$runner" ]]; then
  echo "Missing runner: $runner (ensure it exists and is executable)" >&2
  exit 2
fi

# Discover PARIS IDs
if (( $# >= 3 )); then
  shift 2
  ALIAS_DIR=$1; shift || true
fi

if (( $# > 0 )); then
  PARIS_IDS=("$@")
else
  mapfile -t PARIS_IDS < <(cd "$ALIAS_DIR" && ls -1 energies_list_PARIS*.mac 2>/dev/null | sed -E 's/^energies_list_(PARIS[0-9]+)\.mac$/\1/' | sort -V)
fi

if (( ${#PARIS_IDS[@]} == 0 )); then
  echo "No alias macros found in $ALIAS_DIR" >&2
  exit 0
fi

echo "Found ${#PARIS_IDS[@]} PARIS sets: ${PARIS_IDS[*]}" >&2

for PID in "${PARIS_IDS[@]}"; do
  alias_mac="$ALIAS_DIR/energies_list_${PID}.mac"
  if [[ ! -f "$alias_mac" ]]; then
    echo "Skip $PID (no file $alias_mac)" >&2
    continue
  fi
  echo "=== Running $PID from $alias_mac (EVENTS=$EVENTS, MAX_PROCS=$MAX_PROCS) ===" >&2
  QUIET=${QUIET:-1} G4APP=${G4APP:-./simTetra} "$runner" "$PID" "$alias_mac" "$EVENTS" "$MAX_PROCS"
  echo "=== Done $PID ===" >&2
  echo >&2
done

echo "All requested PARIS runs complete." >&2
