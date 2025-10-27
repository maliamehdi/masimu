#!/bin/bash

# Optimisations possibles pour accélérer les simulations

echo "=== Optimisations suggérées ==="

echo "1. Vérification de l'utilisation disque :"
df -h /srv/data/Geant4/masimu-1/simtetra-main/myanalyse/

echo ""
echo "2. Configuration I/O suggérée :"
echo "   - Utiliser un SSD pour les fichiers ROOT temporaires"
echo "   - Considérer un ramdisk pour les fichiers intermédiaires"

echo ""
echo "3. Monitoring suggéré pendant la simulation :"
echo "   - htop (pour surveiller CPU/RAM)"
echo "   - iotop (pour surveiller I/O disque)"
echo "   - watch 'ls -lah results_*/'"

echo ""
echo "4. Configuration alternative plus agressive :"
echo "   NJOBS=8    # Plus de jobs par angle"
echo "   NTHREADS=5 # Moins de threads par job"  
echo "   MAX_PROCS=8 # Plus de processus concurrents"
echo "   # Total: 8×5=40 threads, mais plus de parallélisme"