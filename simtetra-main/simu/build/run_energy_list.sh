#!/usr/bin/env bash
set -euo pipefail
# Convenience wrapper: take a PARIS alias macro and run one energy per process
# Usage:
#   ./run_energy_list.sh PARIS50 ../bincenters_out_res/energies_list_PARIS50.mac 5000 8
# Args:
#   PARIS_ID    : e.g. PARIS50
#   ALIAS_MAC   : path to energies_list_PARIS50.mac (with /control/alias Elist { ... })
#   EVENTS      : events per energy (/run/beamOn N)
#   MAX_PROCS   : concurrent processes
# Env:
#   G4APP=./simTetra, QUIET=1 to suppress logs

PARIS_ID=${1:?Need PARIS_ID}
ALIAS_MAC=${2:?Need alias macro path}
EVENTS=${3:?Need events per energy}
MAX_PROCS=${4:?Need concurrency}

# Start wall-clock timer
START_TS=$(date +%s)

here=$(cd "$(dirname "$0")" && pwd)
ENERGY_LIST=$(mktemp)
trap 'rm -f "$ENERGY_LIST"' EXIT

# 1) Extract energies to a list
"$here"/extract_alias_to_list.sh "$ALIAS_MAC" "$ENERGY_LIST" >/dev/null

# 2) Run batch (creates per-energy temp macros and deletes them)
G4APP=${G4APP:-./simTetra} QUIET=${QUIET:-1} "$here"/run_single_energy_batch.sh "$PARIS_ID" "$ENERGY_LIST" "$EVENTS" "$MAX_PROCS"

# 3) Done
outdir="../../myanalyse/${PARIS_ID}"
END_TS=$(date +%s)
DUR=$((END_TS-START_TS))
H=$((DUR/3600))
M=$(((DUR%3600)/60))
S=$((DUR%60))
printf "Done. Outputs in %s (one ROOT per energy).\n" "$outdir" >&2
printf "Total wall time: %02d:%02d:%02d (%ds)\n" "$H" "$M" "$S" "$DUR" >&2
