#!/bin/bash

# Script pour renommer et modifier tous les fichiers pompt_PARIS*.mac en prompt_PARIS*.mac
# avec le nouveau contenu template

echo "=== Modification et renommage des fichiers pompt_PARIS*.mac ==="

cd build

# Parcourir tous les fichiers pompt_PARIS*.mac
for OLD_FILE in pompt_PARIS*.mac; do
    if [ -f "$OLD_FILE" ]; then
        # Extraire le numéro PARIS du nom du fichier (ex: pompt_PARIS90.mac -> 90)
        PARIS_NUM=$(echo "$OLD_FILE" | sed 's/pompt_PARIS\([0-9]*\)\.mac/\1/')
        
        # Nouveau nom de fichier
        NEW_FILE="prompt_PARIS${PARIS_NUM}.mac"
        
        echo "→ Traitement de $OLD_FILE -> $NEW_FILE..."
        
        # Créer le nouveau contenu
        cat > "$NEW_FILE" << EOF
# prompt_PARIS${PARIS_NUM}_gun.mac
/process/em/fluo true
/process/em/auger true
/process/em/pixe true

/control/verbose 1
/run/verbose 1
/event/verbose 0
/tracking/verbose 0
/run/initialize

# Particle gun
/gun/particle gamma
/gun/position 0 0 0 mm
/gun/direction 0 0 1

# Paramètres
/control/alias NperE 1000000

# Charge la liste d'énergies { ... } que tu as générée en 1
/control/execute ../bincenters_out/energies_list_PARIS${PARIS_NUM}.mac

# Boucle sur la liste
/control/foreach EkeV \${Elist}
  /gun/energy \${EkeV} keV
  /run/beamOn \${NperE}
/control/endforeach
EOF
        
        # Supprimer l'ancien fichier
        rm "$OLD_FILE"
        
        echo "  ✓ $NEW_FILE créé et $OLD_FILE supprimé"
    fi
done

echo ""
echo "=== Modification terminée ==="
echo "Fichiers prompt_PARIS*.mac créés :"
ls -la prompt_PARIS*.mac 2>/dev/null || echo "Aucun fichier créé"