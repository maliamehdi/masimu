#!/bin/bash

# Script pour créer des fichiers .mac qui lisent directement les énergies
# depuis les fichiers bincenters_PARIS*.txt ligne par ligne

echo "=== Génération des fichiers prompt avec lecture directe ==="

cd build

# Parcourir tous les fichiers bincenters_PARIS*.txt
for BINFILE in ../bincenters_PARIS*.txt; do
    if [ -f "$BINFILE" ]; then
        # Extraire le numéro PARIS du nom du fichier
        PARIS_NUM=$(echo "$BINFILE" | sed 's/.*bincenters_PARIS\([0-9]*\)\.txt/\1/')
        
        # Nom du fichier de sortie
        NEW_FILE="prompt_PARIS${PARIS_NUM}.mac"
        
        echo "→ Création de $NEW_FILE..."
        
        # Créer le header du fichier
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

# Boucle sur toutes les énergies du fichier bincenters
EOF

        # Ajouter une ligne /gun/energy et /run/beamOn pour chaque énergie
        while read -r energy; do
            if [[ "$energy" =~ ^[0-9]+$ ]]; then  # Vérifier que c'est un nombre
                echo "/gun/energy $energy keV" >> "$NEW_FILE"
                echo "/run/beamOn \${NperE}" >> "$NEW_FILE"
            fi
        done < "$BINFILE"
        
        echo "  ✓ $NEW_FILE créé avec $(wc -l < "$BINFILE") énergies"
    fi
done

echo ""
echo "=== Génération terminée ==="
echo "Fichiers créés :"
ls -la prompt_PARIS*.mac 2>/dev/null || echo "Aucun fichier généré"