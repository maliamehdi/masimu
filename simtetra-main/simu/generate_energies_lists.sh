#!/bin/bash

# Script pour générer les fichiers energies_list.mac pour tous les PARIS
# Basé sur le modèle fourni pour PARIS90

echo "=== Génération des listes d'énergies pour tous les PARIS ==="

# Créer le dossier de sortie s'il n'existe pas
mkdir -p bincenters_out

# Parcourir tous les fichiers bincenters_PARIS*.txt
for FILE in bincenters_PARIS*.txt; do
    if [ -f "$FILE" ]; then
        # Extraire le numéro PARIS du nom du fichier (ex: PARIS90 -> 90)
        PARIS_NUM=$(echo "$FILE" | sed 's/bincenters_PARIS\([0-9]*\)\.txt/\1/')
        
        # Nom du fichier de sortie
        OUT="bincenters_out/energies_list_PARIS${PARIS_NUM}.mac"
        
        echo "→ Traitement de $FILE..."
        
        # Générer le fichier .mac
        printf "/control/alias Elist {" > "$OUT"
        tr '\n' ' ' < "$FILE" >> "$OUT"
        printf "}\n" >> "$OUT"
        
        echo "  ✓ Liste macro écrite dans $OUT"
    fi
done

echo ""
echo "=== Génération terminée ==="
echo "Fichiers créés dans le dossier bincenters_out/ :"
ls -la bincenters_out/energies_list_PARIS*.mac 2>/dev/null || echo "Aucun fichier généré"