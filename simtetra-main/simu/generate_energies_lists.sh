#!/usr/bin/env bash
set -euo pipefail

# Script pour générer les fichiers energies_list.mac pour tous les PARIS
# Lit les fichiers bincenters423_out/bincenters_PARIS*.txt et crée un alias
# /control/alias Elist { E1 E2 E3 ... }

echo "=== Génération des listes d'énergies pour tous les PARIS ==="

outdir="bincenters_out_res"
mkdir -p "$outdir"

shopt -s nullglob
files=(bincenters423_out/bincenters_PARIS*.txt)
if (( ${#files[@]} == 0 )); then
    echo "Aucun fichier source trouvé dans bincenters423_out/" >&2
    exit 0
fi

for FILE in "${files[@]}"; do
    if [[ -f "$FILE" ]]; then
        base=$(basename "$FILE")                 # ex: bincenters_PARIS110.txt
        # Extraire le numéro PARIS de manière robuste
        if [[ $base =~ ^bincenters_PARIS([0-9]+)\.txt$ ]]; then
            PARIS_NUM="${BASH_REMATCH[1]}"
        else
            echo "Nom inattendu: $base (ignoré)" >&2
            continue
        fi

        OUT="$outdir/energies_list_PARIS${PARIS_NUM}.mac"
        echo "→ Traitement de $FILE → $OUT"

        # Construire la liste d'énergies sur une seule ligne
        energies=$(paste -sd ' ' "$FILE")
        printf "/control/alias Elist { %s }\n" "$energies" > "$OUT"
        echo "  ✓ Fichier écrit (${OUT})"
    fi
done

echo "\n=== Génération terminée ==="
echo "Fichiers créés dans $outdir :"
ls -la "$outdir"/energies_list_PARIS*.mac 2>/dev/null || echo "Aucun fichier généré"