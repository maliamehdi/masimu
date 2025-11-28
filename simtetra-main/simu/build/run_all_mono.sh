#!/usr/bin/env bash
set -euo pipefail

# Run all mono-energy macros (mono_n_01..mono_n_14) sequentially with
# optimized threads. Organizes outputs into myanalyse/mono by using
# TAG/OUTDIR so Geant4 writes directly to distinct paths per run.
#
# Usage:
#   ./run_all_mono.sh [-j NTHREADS] [--timestamp] [--build DIR]
#
# Options:
#   -j NTHREADS   Number of threads to use (default: nproc)
#   --timestamp   Append YYYYmmdd_HHMMSS to TAG to avoid overwrite
#   --build DIR   Path to build directory containing simTetra
#                 (default: ./build relative to this script)
#
# Notes:
# - Threads are provided via both CLI arg and G4NUM_THREADS to ensure
#   Geant4 MT picks it up early. Macros may set /run/numberOfThreads,
#   but the CLI/env settings are applied earlier and are preferred here.
# - We run from the build directory so relative paths used by the code resolve correctly.
# - We set:
#   * OUTDIR to choose subfolder under myanalyse/
#   * TAG to choose output filename: output_${TAG}.root
#   Paths will look like: myanalyse/${OUTDIR}/output_${TAG}.root

HERE=$(cd "$(dirname "$0")" && pwd)
# Auto-detect layout so the script works from simu/ or simu/build/
if [[ -x "$HERE/simTetra" ]]; then
  # Script placed in build/ next to the binary
  ROOT_DIR="$(cd "$HERE/.." && pwd)"   # points to simu/
  BUILD_DIR="$HERE"
else
  # Script placed in simu/
  ROOT_DIR="$HERE"
  BUILD_DIR="$HERE/build"
fi
NTHREADS="$(command -v nproc >/dev/null 2>&1 && nproc || echo 8)"
USE_TS=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -j)
      NTHREADS="$2"; shift 2;;
    --timestamp)
      USE_TS=1; shift;;
    --build)
      BUILD_DIR="$2"; shift 2;;
    -h|--help)
      sed -n '1,60p' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
    *)
      echo "Unknown option: $1" >&2; exit 2;;
  esac
done

BIN="$BUILD_DIR/simTetra"
if [[ ! -x "$BIN" ]]; then
  echo "Error: simTetra not found at $BIN (build first)" >&2
  exit 1
fi

ts() { date +%Y%m%d_%H%M%S; }
STAMP=""
if [[ $USE_TS -eq 1 ]]; then
  STAMP="_$(ts)"
fi

echo "Using $NTHREADS thread(s). Binary: $BIN"

# Ensure we execute from the build directory so relative paths in the code are valid
pushd "$BUILD_DIR" >/dev/null

# Helper: run a macro by base name using OUTDIR/TAG so RunAction writes to dedicated path
run_with_tag() {
  local base="$1"            # e.g. mono_n_01 (macro basename without .mac)
  local outdir_name="$2"     # e.g. mono (subfolder under myanalyse/)
  local tag_base="$3"        # e.g. mono_n_1 (filename stem)

  local macro_path="$ROOT_DIR/${base}.mac"
  if [[ ! -f "$macro_path" ]]; then
    macro_path="$BUILD_DIR/${base}.mac"
  fi
  if [[ ! -f "$macro_path" ]]; then
    echo "Warning: macro ${base}.mac not found in ROOT_DIR or BUILD_DIR, skipping" >&2
    return 0
  fi

  echo "\n=== Running macro ${base}.mac -> myanalyse/${outdir_name}/output_${tag_base}${STAMP}.root ==="
  local status=0
  (
    export G4NUM_THREADS="$NTHREADS"
    export OUTDIR="$outdir_name"
    export TAG="${tag_base}${STAMP}"
    "$BIN" "$macro_path" "$NTHREADS"
  ) || status=$?
  if [[ $status -ne 0 ]]; then
    echo "Run for macro ${base}.mac failed with code $status" >&2
    return $status
  fi
}

status=0
for n in $(seq -w 1 14); do
  n_short=${n#0}
  base="mono_n_${n}"
  run_with_tag "$base" "mono" "mono_n_${n_short}" || status=$?
  if [[ $status -ne 0 ]]; then
    echo "Run for mono energy index ${n} failed with code $status" >&2
    break
  fi
done

# Also run keV points
if [[ $status -eq 0 ]]; then
  run_with_tag mono_n_100keV mono mono_n_100keV || status=$?
fi
if [[ $status -eq 0 ]]; then
  run_with_tag mono_n_500keV mono mono_n_500keV || status=$?
fi

popd >/dev/null
exit $status
