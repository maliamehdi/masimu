#!/usr/bin/env bash
set -euo pipefail

# Repeat fixed multiplicity runs fixed_mult_01..fixed_mult_11 N times each,
# with different /random/setSeeds each iteration, and progressive hadd.
#
# Usage:
#   ./repeat_fixed_mult_01_to_11.sh 20
#   NTHREADS=36 BASE_SEED=12345 ./repeat_fixed_mult_01_to_11.sh 20
#
# Output:
#   ../../myanalyse/mult<M>/runs/output_fixed_mult_<MM>_it###.root
#   ../../myanalyse/mult<M>/merged.root   (progressive)

NLOOPS=${1:-10}

HERE=$(cd "$(dirname "$0")" && pwd)

# Auto-detect layout so the script works from simu/ or simu/build/
if [[ -x "$HERE/simTetra" ]]; then
  ROOT_DIR="$(cd "$HERE/.." && pwd)"   # simu/
  BUILD_DIR="$HERE"                    # simu/build
else
  ROOT_DIR="$HERE"
  BUILD_DIR="$HERE/build"
fi

BIN="$BUILD_DIR/simTetra"
if [[ ! -x "$BIN" ]]; then
  echo "Error: simTetra not found at $BIN (build first)" >&2
  exit 1
fi

NTHREADS=${NTHREADS:-$(command -v nproc >/dev/null 2>&1 && nproc || echo 8)}
BASE_SEED=${BASE_SEED:-$(date +%s)}

echo "[INFO] BIN=$BIN"
echo "[INFO] ROOT_DIR=$ROOT_DIR"
echo "[INFO] BUILD_DIR=$BUILD_DIR"
echo "[INFO] NTHREADS=$NTHREADS"
echo "[INFO] NLOOPS=$NLOOPS"
echo "[INFO] BASE_SEED=$BASE_SEED"

# Ensure we execute from the build directory so relative paths in the code are valid
pushd "$BUILD_DIR" >/dev/null

find_latest_output_for_base() {
  local base="$1"   # fixed_mult_01 etc
  ls -t ../../myanalyse/output_${base}_run*.root 2>/dev/null | head -n1 || true
}

run_repeat_one_mult() {
  local mm="$1"                # "01".."11"
  local m_int="$2"             # 1..11 (for folder naming & seeding)
  local base="fixed_mult_${mm}"
  local macro="$ROOT_DIR/${base}.mac"

  if [[ ! -f "$macro" ]]; then
    echo "[WARN] missing macro $macro, skipping" >&2
    return 0
  fi

  local dest_dir="../../myanalyse/test_mult${m_int}"
  local runs_dir="${dest_dir}"
  local merged="${dest_dir}/merged.root"
  mkdir -p "$runs_dir"

  echo "===================================================="
  echo "[INFO] Multiplicity ${m_int} (macro: ${base}.mac)"
  echo "[INFO] dest_dir=$dest_dir"

  for ((k=1; k<=NLOOPS; k++)); do
    # Unique deterministic seeds per (m, k)
    s1=$((BASE_SEED + 100000*m_int + 1000*k + 11))
    s2=$((BASE_SEED + 100000*m_int + 1000*k + 97))

    tmp="tmp_${base}_it$(printf "%03d" "$k").mac"
    cat > "$tmp" <<EOF
/control/verbose 0
/run/verbose 0
/event/verbose 0
/tracking/verbose 0
/run/numberOfThreads 36
/random/setSeeds ${s1} ${s2}
/control/execute ${macro}
EOF

    echo "[RUN] ${base}  it=$(printf "%03d" "$k")  seeds=(${s1},${s2})"
    status=0
    (
      export G4NUM_THREADS="$NTHREADS"
      "$BIN" "$tmp" "$NTHREADS"
    ) >/dev/null 2>&1 || status=$?
    rm -f "$tmp"

    if [[ $status -ne 0 ]]; then
      echo "[ERROR] ${base} iteration $k failed (rc=$status)" >&2
      exit $status
    fi

    src=$(find_latest_output_for_base "$base")
    if [[ -z "$src" ]]; then
      echo "[WARN] could not find output for ${base} in ../../myanalyse (pattern output_${base}_run*.root)" >&2
      continue
    fi

    dst="${runs_dir}/output_${base}_it$(printf "%03d" "$k").root"
    mv -f "$src" "$dst"
    echo "[OK] -> $dst"

    # Progressive merge per multiplicity
    if [[ ! -f "$merged" ]]; then
      cp -f "$dst" "$merged"
      echo "[MERGE] init -> $merged"
    else
      tmpm="${dest_dir}/.merged_tmp.root"
      hadd -f "$tmpm" "$merged" "$dst" >/dev/null
      mv -f "$tmpm" "$merged"
      echo "[MERGE] update -> $merged"
    fi
  done

  echo "[DONE] mult${m_int}: merged -> $merged"
}

# Loop multiplicities 01..11
for mm in $(seq -w 11 11); do
  m_int=${mm#0}    # 01->1
  run_repeat_one_mult "$mm" "$m_int"
done

popd >/dev/null
echo "===================================================="
echo "[ALL DONE] fixed_mult_01..11 repeated NLOOPS=$NLOOPS"