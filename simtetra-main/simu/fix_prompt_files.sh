#!/bin/bash

# Script pour créer des fichiers .mac avec boucles d'énergies par plages
# au lieu d'utiliser des alias trop longs

echo "=== Génération des fichiers prompt avec boucles par plages ==="

cd build

# Parcourir tous les fichiers bincenters_PARIS*.txt
for BINFILE in ../bincenters_PARIS*.txt; do
    if [ -f "$BINFILE" ]; then
        # Extraire le numéro PARIS du nom du fichier
        PARIS_NUM=$(echo "$BINFILE" | sed 's/.*bincenters_PARIS\([0-9]*\)\.txt/\1/')
        
        # Nom du fichier de sortie
        NEW_FILE="prompt_PARIS${PARIS_NUM}.mac"
        
        echo "→ Création de $NEW_FILE..."
        
        # Créer le nouveau contenu avec boucle directe
        cat > "$NEW_FILE" << 'EOF'
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

# Boucle directe sur les énergies de 5 à 20000 keV par pas de 10
/control/loop energies_loop.mac EkeV 5 20000 10
EOF

        # Remplacer le placeholder ${PARIS_NUM} par la vraie valeur
        sed -i "s/\${PARIS_NUM}/$PARIS_NUM/g" "$NEW_FILE"
        
        echo "  ✓ $NEW_FILE créé"
    fi
done

# Créer le fichier de boucle énergies_loop.mac
cat > "energies_loop.mac" << 'EOF'
# Fichier de boucle pour les énergies
/gun/energy {EkeV} keV
/run/beamOn {NperE}
EOF

echo ""
echo "=== Génération terminée ==="
echo "Fichiers créés :"
ls -la prompt_PARIS*.mac energies_loop.mac 2>/dev/null || echo "Aucun fichier généré"