#!/usr/bin/env bash
# Run each energy as an independent Geant4 run (one /run/beamOn per process)
# This avoids multi-run file overwrite/corruption issues.
# Usage:
#   ./run_single_energy_batch.sh PARIS50 energies.txt 5000 8
# Args:
#   PARIS_ID      e.g. PARIS50
#   ENERGY_FILE   text file with one energy value in keV per line (e.g. 50, 500, 2000)
#   EVENTS        number of events per energy (beamOn)
#   MAX_PROCS     (optional) concurrent processes (default 1)
# Env overrides:
#   G4APP=./simTetra   path to executable (default ./simTetra)
#   QUIET=1            suppress Geant4 banner
#   NTHREADS=32        Geant4 MT threads per run (default 32)
# Output:
#   ROOT files moved to: ../../myanalyse/<PARIS_ID>/output_<PARIS_ID>_E<energy>keV.root
#   Simple logs: logs/<PARIS_ID>_E<energy>.log (only if QUIET not set)
set -euo pipefail
PARIS_ID=${1:?Need PARIS_ID}
ENERGY_FILE=${2:?Need energy file}
EVENTS=${3:?Need events per energy}
MAX_PROCS=${4:-1}
G4APP=${G4APP:-./simTetra}
NTHREADS=${NTHREADS:-32}
BUILD_DIR=$(pwd)
LOGDIR=logs_single_energy
mkdir -p "$LOGDIR"
ENERGIES=()
# Detect alias-style macro line like: /control/alias Elist {5.5 13.44 ...}
if grep -qE '^[[:space:]]*/control/alias[[:space:]]+[A-Za-z_][A-Za-z0-9_]*[[:space:]]*\{.*\}' "$ENERGY_FILE"; then
  content=$(sed -n 's/.*{\(.*\)}.*/\1/p' "$ENERGY_FILE" | tr -s ' \t' ' ')
  for tok in $content; do
    if [[ $tok =~ ^[0-9]+([.][0-9]+)?$ ]]; then
      ENERGIES+=("$tok")
    fi
  done
else
  while IFS= read -r line; do
    [[ -z "$line" ]] && continue
    # strip comments
    line=${line%%#*}
    line=$(echo "$line" | tr -d ' \t')
    [[ -z "$line" ]] && continue
    # accept pure numbers only
    if [[ $line =~ ^[0-9]+([.][0-9]+)?$ ]]; then
      ENERGIES+=("$line")
    fi
  done < "$ENERGY_FILE"
fi

if (( ${#ENERGIES[@]} == 0 )); then
  echo "No energies parsed from $ENERGY_FILE" >&2
  exit 1
fi

echo "Running ${#ENERGIES[@]} energies for $PARIS_ID with $EVENTS events each (max $MAX_PROCS parallel)" >&2

# Function to run one energy
do_one() {
  local E="$1"; local PARIS_ID="$2"; local EVENTS="$3"; local APP="$4"; local LOGDIR="$5";
  local TAG="${PARIS_ID}_E${E}keV"
  local MAC="tmp_${TAG}.mac"
  local OUTDIR="../../myanalyse/${PARIS_ID}"
  mkdir -p "$OUTDIR"
  # Informative start message
  if [[ -n "${QUIET:-}" ]]; then
    echo "[START][quiet] $TAG" >&2
  else
    echo "[START] $TAG -> $LOGDIR/${TAG}.log" >&2
  fi
  cat > "$MAC" <<EOF
/control/verbose 0
/run/verbose 0
/event/verbose 0
/tracking/verbose 0
/run/numberOfThreads ${NTHREADS}
/run/initialize 
/gun/particle gamma
/gun/position 0 0 -31.8 mm
/gun/energy ${E} keV
/run/beamOn ${EVENTS}
EOF
  if [[ -n "${QUIET:-}" ]]; then
    PARIS_ID="$PARIS_ID" TAG="$TAG" "$APP" "$MAC" >/dev/null 2>&1 || echo "Run failed: $TAG" >&2
  else
    PARIS_ID="$PARIS_ID" TAG="$TAG" "$APP" "$MAC" >"$LOGDIR/${TAG}.log" 2>&1 || echo "Run failed: $TAG" >&2
  fi
  # Move output file into per-PARIS directory if present
  local SRC="../../myanalyse/output_${TAG}.root"
  if [[ -f "$SRC" ]]; then
    mv -f "$SRC" "$OUTDIR/" || true
  fi
  # Done message
  if [[ -n "${QUIET:-}" ]]; then
    echo "[DONE][quiet] $TAG" >&2
  else
    if [[ -f "$OUTDIR/output_${TAG}.root" ]]; then
      echo "[DONE] $TAG -> $OUTDIR/output_${TAG}.root" >&2
    else
      echo "[DONE] $TAG (no output found)" >&2
    fi
  fi
  rm -f "$MAC"
}

# Parallel dispatcher (portable, avoids 'wait -n' requirement)
for E in "${ENERGIES[@]}"; do
  # Throttle to MAX_PROCS concurrent background jobs
  while true; do
    running_jobs=$(jobs -r -p | wc -l)
    if (( running_jobs < MAX_PROCS )); then
      break
    fi
    sleep 0.2
  done
  do_one "$E" "$PARIS_ID" "$EVENTS" "$G4APP" "$LOGDIR" &
done
wait

echo "All runs complete. Output files in ../../myanalyse/${PARIS_ID}/ output_${PARIS_ID}_E*.root" >&2
