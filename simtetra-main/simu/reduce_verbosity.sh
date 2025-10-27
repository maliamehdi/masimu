#!/bin/bash

# Script pour réduire la verbosité des fichiers prompt_PARIS*.mac
# et minimiser les messages dans le terminal

echo "=== Réduction de la verbosité des fichiers prompt_PARIS*.mac ==="

cd build

# Parcourir tous les fichiers prompt_PARIS*.mac
for FILE in prompt_PARIS*.mac; do
    if [ -f "$FILE" ]; then
        echo "→ Modification de $FILE..."
        
        # Créer une sauvegarde
        cp "$FILE" "${FILE}.backup"
        
        # Remplacer les paramètres de verbosité pour minimiser les messages
        sed -i 's|/control/verbose 1|/control/verbose 0|g' "$FILE"
        sed -i 's|/run/verbose 1|/run/verbose 0|g' "$FILE"
        
        # S'assurer que event et tracking verbosity restent à 0
        sed -i 's|/event/verbose.*|/event/verbose 0|g' "$FILE"
        sed -i 's|/tracking/verbose.*|/tracking/verbose 0|g' "$FILE"
        
        echo "  ✓ $FILE modifié (verbosité réduite)"
    fi
done

echo ""
echo "=== Modification terminée ==="
echo "Nouveaux paramètres de verbosité :"
echo "  /control/verbose 0   (pas de messages de contrôle)"
echo "  /run/verbose 0       (pas de messages de run)"
echo "  /event/verbose 0     (pas de messages d'événements)"
echo "  /tracking/verbose 0  (pas de messages de tracking)"
echo ""
echo "Les sauvegardes sont dans les fichiers *.backup"