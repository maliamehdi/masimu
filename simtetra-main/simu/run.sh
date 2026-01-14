#!/usr/bin/env bash
set -euo pipefail

########################
# Paramètres utilisateur
########################
EXE="./simTetra"                   # binaire
COMMON_MAC="152Eu.mac"      # macro appelée après /run/initialize (ta config source etc.)
NEVT=100000                        # évts par job (2.5M pour finir rapidement)
NTHREADS=8                         # threads Geant4 par job 
NJOBS=10000                           # nombre total de jobs (10 × 2.5M = 25M événements)
MAX_PROCS=1                        # nb maximum de processus concurrents
OUTBASE="eu152bis"                    # préfixe pour nommer les fichiers
RESULTS_DIR="../../myanalyse/results_$(date +%Y%m%d_%H%M%S)"  # dossier où ranger logs et ROOT fusionnés
# Dossier où ton code écrit les ROOT (chemin utilisé par MyRunAction)
ANA_DIR="../../myanalyse"

# Graine de base (unique par lot) pour dériver des seeds stables
BASE_SEED=$(date +%s)

mkdir -p "${RESULTS_DIR}"
echo "[INFO] Résultats dans: ${RESULTS_DIR}"

# Fonction qui fabrique la macro de run "run_i.mac"
make_macro() {
  local i=$1
  local mac="run_${i}.mac"
  local s1=$((BASE_SEED + 1000*i + 123))
  local s2=$((BASE_SEED + 1000*i + 987))

  cat > "${mac}" <<EOF
/run/numberOfThreads ${NTHREADS}
/random/setSeeds ${s1} ${s2}

# -- tes options RDM/EM, etc. si besoin --
#/process/em/fluo true
#/process/em/auger true
#/process/em/pixe true
#/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+60 year
#/gun/particle ion
#/gun/ion 63 152 0 0
#/gun/energy 0 keV
#/gun/position 0 0 0 mm

/run/initialize
/control/execute ${COMMON_MAC}
/run/beamOn ${NEVT}
EOF
  echo "${mac}"
}

# Fonction qui lance 1 job + déplace son ROOT de ../../myanalyse vers RESULTS_DIR
run_one() {
  local i=$1
  local tag="${OUTBASE}_r$(printf "%03d" ${i})"
  local mac
  mac=$(make_macro "${i}")

  echo "[JOB ${i}] start : TAG=${tag}, macro=${mac}"
  # Export TAG pour que MyRunAction nomme le fichier root en conséquence
  TAG="${tag}" "${EXE}" "${mac}" > "${RESULTS_DIR}/log_${i}.txt" 2>&1

  # Le code écrit normalement dans ${ANA_DIR}/output_${TAG}.root
  local src="${ANA_DIR}/output_${tag}.root"
  local dst="${RESULTS_DIR}/${OUTBASE}_$(printf "%03d" ${i}).root"

  if [[ -f "${src}" ]]; then
    mv -f "${src}" "${dst}"
    echo "[JOB ${i}] OK : ${dst}"
  else
    echo "[JOB ${i}] WARN : fichier ${src} introuvable (vérifie le chemin dans MyRunAction)."
  fi
}

# Limiteur simple de parallélisme (<= MAX_PROCS processus à la fois)
wait_for_slot() {
  while (( $(jobs -rp | wc -l) >= MAX_PROCS )); do
    sleep 0.2
  done
}

########################
# Lancement des jobs
########################
echo "[INFO] NJOBS=${NJOBS}, NEVT=${NEVT}/job, NTHREADS=${NTHREADS}/job, MAX_PROCS=${MAX_PROCS}"

for i in $(seq 0 $((NJOBS-1))); do
  wait_for_slot
  run_one "${i}" &
done

wait
echo "[INFO] Tous les jobs sont terminés."

########################
# Fusion avec hadd
########################
MERGED="${RESULTS_DIR}/${OUTBASE}_merged.root"
echo "[INFO] Fusion: ${MERGED}"
# -f : force overwrite si le fichier existe
hadd -f "${MERGED}" "${RESULTS_DIR}"/${OUTBASE}_*.root

echo "[DONE] Fichier fusionné: ${MERGED}"