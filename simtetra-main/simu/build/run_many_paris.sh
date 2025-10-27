#!/usr/bin/env bash
set -euo pipefail

########################
# Paramètres utilisateur
########################
EXE="./simTetra"                         # binaire Geant4
PARIS_IDS=(50 70 90 110 130 235 262 278 305)  # angles
NJOBS=4                                  # nb de jobs PAR angle
NTHREADS=10                              # threads Geant4 par job (optimisé pour 40 cœurs)
MAX_PROCS=4                              # nb max de processus concurrents (4×10=40 threads)
OUTBASE="paris"                          # préfixe pour les tags ROOT
RESULTS_DIR="../../myanalyse/results_$(date +%Y%m%d_%H%M%S)"
ANA_DIR="../../myanalyse"                # dossier où MyRunAction écrit les ROOT

# graine de base pour rendre les seeds reproductibles par lot
BASE_SEED=$(date +%s)

mkdir -p "${RESULTS_DIR}"
echo "[INFO] Résultats dans: ${RESULTS_DIR}"

make_wrapper_macro() {
  local id="$1"     # angle PARIS (ex: 70)
  local i="$2"      # index de job (0..NJOBS-1)
  local mac="run_${id}_${i}.mac"

  # deux seeds dérivées du lot / angle / job
  local s1=$((BASE_SEED + 100000*id + 1000*i + 123))
  local s2=$((BASE_SEED + 100000*id + 1000*i + 987))

  cat > "${mac}" <<EOF
/run/numberOfThreads ${NTHREADS}
/random/setSeeds ${s1} ${s2}

# IMPORTANT: on configure avant l'initialize de la macro PARIS
/control/execute prompt_PARIS${id}.mac
EOF
  echo "${mac}"
}

run_one() {
  local id="$1"
  local i="$2"
  local tag="${OUTBASE}_PARIS${id}_r$(printf "%03d" ${i})"

  local mac
  mac=$(make_wrapper_macro "${id}" "${i}")
  echo "[JOB PARIS${id} #${i}] start : TAG=${tag}, macro=${mac}"

  # -> ton RunAction doit utiliser \$TAG pour nommer le fichier ROOT
  TAG="${tag}" "${EXE}" "${mac}" > "${RESULTS_DIR}/log_PARIS${id}_r${i}.txt" 2>&1

  # Le code écrit: ${ANA_DIR}/output_${tag}.root
  local src="${ANA_DIR}/output_${tag}.root"
  local dst="${RESULTS_DIR}/${tag}.root"
  if [[ -f "${src}" ]]; then
    mv -f "${src}" "${dst}"
    echo "[JOB PARIS${id} #${i}] OK : ${dst}"
  else
    echo "[JOB PARIS${id} #${i}] WARN : ${src} introuvable. (Vérifie BeginOfRunAction et TAG)"
  fi
}

wait_for_slot() {
  while (( $(jobs -rp | wc -l) >= MAX_PROCS )); do
    sleep 0.2
  done
}

echo "[INFO] NTHREADS/job=${NTHREADS}, MAX_PROCS=${MAX_PROCS}"
echo "[INFO] Contrainte conseillée: MAX_PROCS * NTHREADS ≤ nb coeurs"

###############
# Lancements
###############
for id in "${PARIS_IDS[@]}"; do
  for i in $(seq 0 $((NJOBS-1))); do
    wait_for_slot
    run_one "${id}" "${i}" &
  done
done

wait
echo "[INFO] Tous les jobs sont terminés."

########################
# Fusion par angle
########################
for id in "${PARIS_IDS[@]}"; do
  pat="${RESULTS_DIR}/${OUTBASE}_PARIS${id}_r"*.root
  if ls ${pat} >/dev/null 2>&1; then
    out="${RESULTS_DIR}/${OUTBASE}_PARIS${id}_merged.root"
    echo "[HADD] PARIS${id} -> ${out}"
    hadd -f "${out}" ${pat}
  else
    echo "[HADD] PARIS${id}: aucun fichier à fusionner"
  fi
done

########################
# Fusion globale (tous angles)
########################
glob_list=()
for id in "${PARIS_IDS[@]}"; do
  f="${RESULTS_DIR}/${OUTBASE}_PARIS${id}_merged.root"
  [[ -f "${f}" ]] && glob_list+=("${f}")
done

if (( ${#glob_list[@]} > 0 )); then
  MERGED_ALL="${RESULTS_DIR}/${OUTBASE}_ALLangles_merged.root"
  echo "[HADD] global -> ${MERGED_ALL}"
  hadd -f "${MERGED_ALL}" "${glob_list[@]}"
  echo "[DONE] Fusion globale: ${MERGED_ALL}"
else
  echo "[HADD] global: rien à fusionner"
fi
# Optionnel : supprimer les ROOT par job
rm -f "${RESULTS_DIR}"/paris_PARIS*_r*.root