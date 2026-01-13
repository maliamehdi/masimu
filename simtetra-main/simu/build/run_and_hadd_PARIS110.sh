#!/usr/bin/env bash
set -euo pipefail

SIMU_DIR="/srv/data/Geant4/masimu-1/simtetra-main/simu/build"
ANA_DIR="/srv/data/Geant4/masimu-1/simtetra-main/myanalyse/PARIS110"
RUN_CMD="./run_paris.scr"

FILES_TO_MERGE=(
  "output_PARIS110_E64.7369keV.root"
  "output_PARIS110_E1741.87keV.root"
  "output_PARIS110_E7554.07keV.root"
  "output_PARIS110_E11734.5keV.root"
)

NLOOPS="${1:-20}"

die() { echo "[ERROR] $*" >&2; exit 1; }
check_file() { [[ -f "$1" ]] || die "Missing file: $ANA_DIR/$1"; }

echo "[INFO] Will run ${NLOOPS} iterations."

# ------------------------------------------------------------
# If complement1.root doesn't exist, create it from the 4 files
# ------------------------------------------------------------
cd "$ANA_DIR"
if [[ ! -f "complement1.root" ]]; then
  echo "[INFO] complement1.root not found -> will create it after first ./run_paris.scr"
fi

for ((k=1; k<=NLOOPS; k++)); do
  echo "===================================================="
  echo "[INFO] Iteration ${k}/${NLOOPS}"

  # 1) Run simulation
  cd "$SIMU_DIR"
  echo "[INFO] Running: $RUN_CMD"
  $RUN_CMD

  # 2) Merge outputs
  cd "$ANA_DIR"

  # check regenerated files
  for f in "${FILES_TO_MERGE[@]}"; do
    check_file "$f"
  done

  if [[ ! -f "complement1.root" ]]; then
    # First complement: hadd of the 4 files only
    echo "[INFO] Creating complement1.root from the 4 output files"
    hadd -f "complement1.root" "${FILES_TO_MERGE[@]}"
    echo "[INFO] Done -> $ANA_DIR/complement1.root"
    continue
  fi

  # If complement1 exists, we append iteration-wise:
  # complement(k+1) = hadd(complement(k), 4 files)
  prev_idx=$k
  next_idx=$((k+1))

  prev="complement${prev_idx}.root"
  out="complement${next_idx}.root"

  # If prev doesn't exist (e.g. you already had complement1 but not complement2),
  # we start from the last available complement instead.
  if [[ ! -f "$prev" ]]; then
    # find last existing complementN.root
    last=$(ls -1 complement*.root 2>/dev/null | sed -E 's/^complement([0-9]+)\.root$/\1/' | sort -n | tail -1 || true)
    [[ -n "${last:-}" ]] || die "No complement*.root found but complement1.root check passed?"
    prev="complement${last}.root"
    out="complement$((last+1)).root"
  fi

  [[ -f "$prev" ]] || die "Missing previous complement: $prev"
  [[ ! -f "$out" ]] || die "Output already exists: $out (refusing to overwrite)"

  echo "[INFO] Creating $out from:"
  echo "       - $prev"
  printf "       - %s\n" "${FILES_TO_MERGE[@]}"

  hadd -f "$out" "$prev" "${FILES_TO_MERGE[@]}"

  echo "[INFO] Done -> $ANA_DIR/$out"
done

echo "===================================================="
echo "[INFO] All done."