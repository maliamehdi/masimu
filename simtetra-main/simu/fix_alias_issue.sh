#!/bin/bash

# Script pour corriger le problème des alias dans les fichiers prompt_PARIS*.mac
# Remplace ${NperE} par la valeur directe 1000000

echo "=== Correction des fichiers prompt_PARIS*.mac ==="

cd build

# Parcourir tous les fichiers prompt_PARIS*.mac
for FILE in prompt_PARIS*.mac; do
    if [ -f "$FILE" ]; then
        echo "→ Correction de $FILE..."
        
        # Remplacer ${NperE} par 1000000 directement
        sed -i 's/\${NperE}/1000000/g' "$FILE"
        
        # Supprimer la ligne d'alias qui n'est plus nécessaire
        sed -i '/\/control\/alias NperE 1000000/d' "$FILE"
        
        echo "  ✓ $FILE corrigé"
    fi
done

echo ""
echo "=== Correction terminée ==="
echo "Les fichiers ont été corrigés :"
ls -la prompt_PARIS*.mac 2>/dev/null || echo "Aucun fichier trouvé"