#!/usr/bin/env bash
set -euo pipefail

########################
# Paramètres utilisateur
########################
EXE="./simTetra"                              # binaire Geant4
PARIS_IDS=(50 70 90 110 130 235 262 278 305) # angles
NJOBS=${NJOBS:-100}                           # nb de jobs par angle
NTHREADS=${NTHREADS:-10}                      # threads Geant4 par job
MAX_PROCS=${MAX_PROCS:-2}                     # nb max de processus concurrents
OUTBASE="paris"                               # préfixe des tags ROOT
RESULTS_DIR="../../myanalyse/results_$(date +%Y%m%d_%H%M%S)"
ANA_DIR="../../myanalyse"                     # dossier où MyRunAction écrit les ROOT
NREPEAT=${NREPEAT:-1}                         # nb de répétitions d’une même macro dans un job

# Logs: SAVE_LOGS=0 pour /dev/null, sinon fichier; GZIP_LOGS=1 pour compresser les logs
SAVE_LOGS=${SAVE_LOGS:-1}
GZIP_LOGS=${GZIP_LOGS:-0}

# Auto-ajustement optionnel de la concurrence (export AUTO_CONCURRENCY=1 pour activer)
AUTO_CONCURRENCY=${AUTO_CONCURRENCY:-1}
if [[ "${AUTO_CONCURRENCY}" == "1" ]]; then
  RESERVE_CORES=${RESERVE_CORES:-4}
  TOTAL_CORES=$(nproc)
  AVAIL=$(( TOTAL_CORES - RESERVE_CORES ))
  (( AVAIL < 1 )) && AVAIL=1
  TMP=$(( AVAIL / NTHREADS ))
  (( TMP < 1 )) && TMP=1
  MAX_PROCS=${TMP}
fi

# graine de base pour rendre les seeds reproductibles par lot
BASE_SEED=$(date +%s)

mkdir -p "${RESULTS_DIR}"
echo "[INFO] Résultats dans: ${RESULTS_DIR}"
echo "[INFO] NTHREADS/job=${NTHREADS}, MAX_PROCS=${MAX_PROCS}, NREPEAT/job=${NREPEAT}"
echo "[INFO] Contrainte conseillée: MAX_PROCS * NTHREADS ≤ nproc=$(nproc)"

########################################
# Génération d'une macro wrapper par run
########################################
make_wrapper_macro() {
  local id="$1"                 # angle PARIS (ex: 70)
  local i="$2"                  # index de job (0..NJOBS-1)
  local rep_offset="${3:-0}"    # répétition (0..NREPEAT-1) pour varier seeds/nom
  local mac="run_${id}_${i}_rep${rep_offset}.mac"

  # deux seeds dérivées du lot / angle / job / répétition
  local s1=$((BASE_SEED + 100000*id + 1000*i + 10*rep_offset + 123))
  local s2=$((BASE_SEED + 100000*id + 1000*i + 10*rep_offset + 987))

  cat > "${mac}" <<EOF
/control/verbose 0
/run/verbose 0
/event/verbose 0
/tracking/verbose 0
/run/numberOfThreads ${NTHREADS}
/random/setSeeds ${s1} ${s2}

# IMPORTANT: configure le setup PARIS avant l'initialize
/control/execute prompt_PARIS${id}.mac
EOF
  echo "${mac}"
}

#######################
# Lancement d'un "job"
#######################
run_one() {
  local id="$1"
  local i="$2"
  local tag_base="${OUTBASE}_PARIS${id}_r$(printf "%03d" ${i})"

  echo "[JOB PARIS${id} #${i}] start : TAG base=${tag_base}, NREPEAT=${NREPEAT}"
  local logfile="${RESULTS_DIR}/log_PARIS${id}_r${i}.txt"

  for rep in $(seq 0 $((NREPEAT-1))); do
    local repTag="${tag_base}_p$(printf "%03d" ${rep})"
    local mac
    mac=$(make_wrapper_macro "${id}" "${i}" "${rep}")

    # Exécute simTetra (RunAction ouvre 1 seul ROOT par exécution grâce à TAG)
    if [[ "${SAVE_LOGS}" == "0" ]]; then
      QUIET=1 TAG="${repTag}" "${EXE}" "${mac}" > /dev/null 2>&1
    else
      QUIET=1 TAG="${repTag}" "${EXE}" "${mac}" >> "${logfile}" 2>&1
    fi

    # Le code écrit: ${ANA_DIR}/output_${repTag}.root (cf. RunAction::BeginOfRunAction)
    local src="${ANA_DIR}/output_${repTag}.root"
    local dst="${RESULTS_DIR}/${repTag}.root"

    if [[ -f "${src}" ]]; then
      mv -f "${src}" "${dst}"
      echo "[JOB PARIS${id} #${i}] saved: ${dst}"
    else
      echo "[JOB PARIS${id} #${i}] WARN : ${src} introuvable."
    fi
  done

  echo "[JOB PARIS${id} #${i}] DONE"
  if [[ "${SAVE_LOGS}" != "0" && "${GZIP_LOGS}" == "1" && -f "${logfile}" ]]; then
    gzip -f "${logfile}" || true
  fi
}

##############################
# Régulation de la concurrence
##############################
wait_for_slot() {
  while (( $(jobs -rp | wc -l) >= MAX_PROCS )); do
    sleep 0.2
  done
}

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
# Fusion par angle (FIN)
########################
for id in "${PARIS_IDS[@]}"; do
  # Tous les ROOT produits pour cet angle (tous jobs + répétitions)
  pat="${RESULTS_DIR}/${OUTBASE}_PARIS${id}_r"*.root
  if ls ${pat} >/dev/null 2>&1; then
    out="${RESULTS_DIR}/${OUTBASE}_PARIS${id}_merged.root"

    # Sauvegarde d'un merged existant s'il existe
    if [[ -f "${out}" ]]; then
      ts=$(date +%Y%m%d_%H%M%S)
      prev="${out%.root}_prev_${ts}.root"
      echo "[SAFE] ${out} existe, sauvegarde -> ${prev}"
      mv -f "${out}" "${prev}"
    fi

    # Log des sources utilisées pour transparence/reproductibilité
    src_list_file="${RESULTS_DIR}/${OUTBASE}_PARIS${id}_sources_$(date +%Y%m%d_%H%M%S).lst"
    ls -1 ${pat} > "${src_list_file}"

    echo "[HADD] PARIS${id} -> ${out} (sources -> ${src_list_file})"
    hadd -f "${out}" ${pat}
  else
    echo "[HADD] PARIS${id}: aucun fichier à fusionner"
  fi
done

echo "[DONE] Fusions par angle terminées."