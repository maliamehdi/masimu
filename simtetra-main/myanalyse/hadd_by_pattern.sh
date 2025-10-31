#!/usr/bin/env bash
set -euo pipefail

# Merge ROOT files whose filenames match a pattern using hadd.
#
# Usage:
#   ./hadd_by_pattern.sh -o output.root -p "pattern" [-d dir] [-r] [--delete] [--dry-run] [--] [extra hadd args]
#
# Examples:
#   ./hadd_by_pattern.sh -o merged.root -p "run42" -d ../myanalyse
#   ./hadd_by_pattern.sh -o cf252_all.root -p "cf252" -r -d ../myanalyse
#   ./hadd_by_pattern.sh -o out.root -p "*.root" -d ../myanalyse -- -k  # pass -k to hadd
#   ./hadd_by_pattern.sh -o merged.root -p "tag" -d ./data --delete     # delete inputs after merging
#
# Notes:
# - Pattern is matched on the basename; by default we treat it as a shell glob if it contains *?[.
#   Otherwise we do a substring match ("*pattern*"). Quote the pattern to avoid shell expansion.
# - Use -r for recursive search under -d. Without -d, current directory is used.
# - Extra args after "--" are passed to hadd (e.g., -k, -f).

out=""
dir="."
pattern=""
recursive=0
dry=0
delete_inputs=0

print_help() {
  sed -n '1,50p' "$0" | sed 's/^# \{0,1\}//'
}

# Parse options
while [[ $# -gt 0 ]]; do
  case "$1" in
    -o|--output)
      out="$2"; shift 2;;
    -d|--dir)
      dir="$2"; shift 2;;
    -p|--pattern)
      pattern="$2"; shift 2;;
    -r|--recursive)
      recursive=1; shift;;
    -x|--delete|--rm)
      delete_inputs=1; shift;;
    --dry-run)
      dry=1; shift;;
    -h|--help)
      print_help; exit 0;;
    --)
      shift; break;;
    *)
      echo "Unknown option: $1" >&2; print_help; exit 2;;
  esac
done

# Remaining args are passed to hadd
HADD_EXTRA_ARGS=("$@")

if [[ -z "$out" || -z "$pattern" ]]; then
  echo "Error: -o <output.root> and -p <pattern> are required" >&2
  print_help; exit 2
fi

if [[ -e "$out" ]]; then
  echo "Warning: output file '$out' exists and will be overwritten by hadd -f (if specified)." >&2
fi

# Build find command
if [[ $recursive -eq 1 ]]; then
  find_cmd=(find "$dir" -type f -name "*.root")
else
  find_cmd=(find "$dir" -maxdepth 1 -type f -name "*.root")
fi

# Decide matching mode: glob-like vs substring
match_glob=0
if [[ "$pattern" == *'*'* || "$pattern" == *'?'* || "$pattern" == *'['* ]]; then
  match_glob=1
fi

mapfile -t files < <("${find_cmd[@]}" -print0 | xargs -0 -n1 basename | awk -v p="$pattern" -v g="$match_glob" '{
  if (g==1) {
    # Faux glob via bash pattern-matching; we'll test later in bash
    print $0
  } else {
    # substring
    if (index($0, p)>0) print $0
  }
}' | sort -u)

# Reconstruct paths
selected=()
out_base=$(basename -- "$out")
if [[ ${#files[@]} -gt 0 ]]; then
  for f in "${files[@]}"; do
    # Apply glob filter in bash if requested
    if [[ $match_glob -eq 1 ]]; then
      if [[ ! $f == $pattern ]]; then continue; fi
    fi
    # Never include the output file itself as input
    if [[ "$f" == "$out_base" ]]; then continue; fi
    # find path in dir (prefer top-level when not recursive)
    if [[ $recursive -eq 1 ]]; then
      # search first match
      pth=$(find "$dir" -type f -name "$f" -print -quit)
    else
      pth="$dir/$f"
    fi
    [[ -f "$pth" ]] && selected+=("$pth")
  done
fi

if [[ ${#selected[@]} -eq 0 ]]; then
  echo "No files found in '$dir' matching pattern '$pattern'" >&2
  exit 1
fi

# Show what we will do
echo "Found ${#selected[@]} file(s):"
for f in "${selected[@]}"; do echo "  $f"; done

# Run hadd
cmd=(hadd "${HADD_EXTRA_ARGS[@]}" "$out")
cmd+=("${selected[@]}")

echo "\nRunning: ${cmd[*]}"
if [[ $dry -eq 1 ]]; then
  if [[ $delete_inputs -eq 1 ]]; then
    echo "(Dry-run) Would delete input files after successful merge:"
    for f in "${selected[@]}"; do echo "  rm -- \"$f\""; done
  fi
  echo "Dry-run only. Exiting."; exit 0
fi

"${cmd[@]}"

echo "Done. Wrote: $out"

# Optionally delete inputs after successful merge
if [[ $delete_inputs -eq 1 ]]; then
  echo "Deleting input files (--delete):"
  for f in "${selected[@]}"; do
    # Double guard: don't delete the output file even if it somehow appears
    if [[ $(basename -- "$f") == "$out_base" ]]; then continue; fi
    rm -f -- "$f"
    echo "  removed: $f"
  done
fi