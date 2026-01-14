#!/usr/bin/env bash

echo "=== LANCEMENT SÉCURISÉ DE LA SIMULATION ==="
echo "Cette simulation résistera à :"
echo "- Verrouillage d'écran ✅"
echo "- Fermeture du terminal ✅" 
echo "- Déconnexion SSH ✅"
echo ""

# Désactiver la mise en veille (nécessite sudo)
echo "1. Désactivation de la mise en veille..."
sudo pmset -c sleep 0 -c displaysleep 60
echo "   ✅ Mise en veille désactivée"

# Lancer avec nohup pour résister à la fermeture de terminal
echo "2. Lancement de la simulation en arrière-plan..."
nohup ./run.sh > simulation.log 2>&1 &
SIMUL_PID=$!

echo "   ✅ Simulation lancée (PID: $SIMUL_PID)"
echo "   📄 Logs dans: simulation.log"

# Créer un fichier de suivi
echo "$SIMUL_PID" > simulation.pid
echo "$(date): Simulation démarrée (PID: $SIMUL_PID)" >> simulation_status.log

echo ""
echo "=== COMMANDES UTILES ==="
echo "• Suivre les logs:     tail -f simulation.log"
echo "• Vérifier l'état:     ps aux | grep simTetra"
echo "• Arrêter si besoin:   kill $SIMUL_PID"
echo "• Progression:         ls -la results_*/"
echo ""
echo "🚀 Vous pouvez maintenant verrouiller votre Mac en toute sécurité !"
echo "   La simulation continuera en arrière-plan."

# Option : surveiller le début
echo "Vérification du démarrage (10 secondes)..."
sleep 10
if ps -p $SIMUL_PID > /dev/null; then
    echo "✅ Simulation en cours d'exécution"
    echo "📊 $(ps -p $SIMUL_PID -o pid,pcpu,pmem,time)"
else
    echo "❌ Problème au démarrage, vérifiez simulation.log"
fi