#!/usr/bin/env bash
set -euo pipefail

# Run all fixed-multiplicity macros (1..12) sequentially with
# Watt energy sampling and optimized threads, then run extra macros:
#   - auto252Cf_neutrons.mac
#   - mono_n.mac
#   - duo.mac
# Organizes outputs into myanalyse/mult1..mult12 without touching C++ code.
#
# Usage:
#   ./run_all_fixed_mult.sh [-j NTHREADS] [--timestamp] [--build DIR]
#
# Options:
#   -j NTHREADS   Number of threads to use (default: nproc)
#   --timestamp   Append YYYYmmdd_HHMMSS to TAG to avoid overwrite
#   --build DIR   Path to build directory containing simTetra
#                 (default: ./build relative to this script)
#
# Notes:
# - Threads are provided via both CLI arg and G4NUM_THREADS to ensure
#   Geant4 MT picks it up early. Macros also set /run/numberOfThreads 8,
#   but the CLI/env settings are applied earlier and are preferred here.
# - We run from the build directory so relative paths used by the code resolve correctly.
# - After each run, we move the produced ROOT file into a tidy subfolder:
#   * fixed_mult_N  -> myanalyse/multN/output_fixed_mult_N[+timestamp].root
#   * auto252Cf     -> myanalyse/auto252Cf/output_auto252Cf[+timestamp].root
#   * mono_n        -> myanalyse/mono/output_mono_n[+timestamp].root
#   * duo           -> myanalyse/duo/output_duo[+timestamp].root

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

# Helper: run a macro by base name and move the latest output into a named folder/file
run_and_move() {
  local base="$1"       # e.g. auto252Cf_neutrons
  local dest_dir_name="$2"  # e.g. auto252Cf
  local dest_basename="$3"  # e.g. auto252Cf

  local macro_path="$ROOT_DIR/${base}.mac"
  if [[ ! -f "$macro_path" ]]; then
    macro_path="$BUILD_DIR/${base}.mac"
  fi
  if [[ ! -f "$macro_path" ]]; then
    echo "Warning: macro ${base}.mac not found in ROOT_DIR or BUILD_DIR, skipping" >&2
    return 0
  fi

  echo "\n=== Running macro ${base}.mac ==="
  local status=0
  (
    export G4NUM_THREADS="$NTHREADS"
    "$BIN" "$macro_path" "$NTHREADS"
  ) || status=$?
  if [[ $status -ne 0 ]]; then
    echo "Run for macro ${base}.mac failed with code $status" >&2
    return $status
  fi

  # Move latest output for this base
  local src
  src=$(ls -t ../../myanalyse/output_${base}_run*.root 2>/dev/null | head -n1 || true)
  if [[ -z "$src" ]]; then
    echo "Warning: could not find output for ${base} in ../../myanalyse" >&2
    return 0
  fi
  local dest_dir="../../myanalyse/${dest_dir_name}"
  mkdir -p "$dest_dir"
  local dest_file="$dest_dir/output_${dest_basename}${STAMP}.root"
  echo "Moving $src -> $dest_file"
  mv -f "$src" "$dest_file"
}

status=0
for n in $(seq -w 11 11); do
  macro="$ROOT_DIR/fixed_mult_${n}.mac"
  if [[ ! -f "$macro" ]]; then
    echo "Warning: missing macro $macro, skipping" >&2
    continue
  fi
  echo "\n=== Running multiplicity ${n} ==="
  (
    export G4NUM_THREADS="$NTHREADS"
    # Pass threads as 2nd CLI arg too (applied before macro is read)
    "$BIN" "$macro" "$NTHREADS"
  ) || status=$?
  if [[ $status -ne 0 ]]; then
    echo "Run for multiplicity ${n} failed with code $status" >&2
    break
  fi

  # After successful run, move the newest produced ROOT file for this macro
  base="fixed_mult_${n}"
  src=$(ls -t ../../myanalyse/output_${base}_run*.root 2>/dev/null | head -n1 || true)
  if [[ -z "$src" ]]; then
    echo "Warning: could not find output for ${base} in ../../myanalyse" >&2
    continue
  fi
  # Convert 01 -> 1 for directory naming
  n_short=${n#0}
  dest_dir=../../myanalyse/mult${n_short}
  mkdir -p "$dest_dir"
  dest_file="$dest_dir/output_fixed_mult_${n_short}${STAMP}.root"
  echo "Moving $src -> $dest_file"
  mv -f "$src" "$dest_file"
done

# Extra macros after the fixed multiplicities
#run_and_move auto252Cf_neutrons auto252Cf auto252Cf || status=$?
#if [[ $status -ne 0 ]]; then popd >/dev/null; exit $status; fi

#run_and_move mono_n mono mono_n || status=$?
#if [[ $status -ne 0 ]]; then popd >/dev/null; exit $status; fi

#run_and_move duo duo duo || status=$?
#if [[ $status -ne 0 ]]; then popd >/dev/null; exit $status; fi

popd >/dev/null
exit $status
