#include <iostream>

#include "G4RunManagerFactory.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4UImanager.hh"
#include "G4SteppingVerbose.hh"
#include "Randomize.hh"
#include <chrono>
#include <ctime>

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"

#include "G4ParticleHPManager.hh"

// Si tu fermes l’analyse ici, dé-commente les deux includes suivants
// #include "g4analysis.hh"     // fabrique d’AnalysisManager (Geant4 11.x)
// #include "G4AutoDelete.hh"

int main(int argc, char** argv)
{
  // --- Run manager (choisit tout seul Serial/MT selon la build)
  auto* runManager = G4RunManagerFactory::CreateRunManager();
  #ifdef G4MULTITHREADED
     runManager->SetNumberOfThreads(G4Threading::G4GetNumberOfCores());
  #endif

  // RNG + verbosité unités
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  G4SteppingVerbose::UseBestUnit(4);

  // Detector / Physics / Actions
  runManager->SetUserInitialization(new MyDetectorConstruction());
  runManager->SetUserInitialization(new MyPhysicsList());
  // Pass the actual macro name (argv[1]) to actions for proper output naming
  G4String macroName = (argc > 1 ? G4String(argv[1]) : G4String("vis.mac"));
  runManager->SetUserInitialization(new MyActionInitialization(macroName));

  // Réglages HP (ok ici, avant /run/initialize)
  auto* hp = G4ParticleHPManager::GetInstance();
  hp->SetSkipMissingIsotopes(true);
  hp->SetDoNotAdjustFinalState(true);
  hp->SetUseOnlyPhotoEvaporation(false);
  hp->SetNeglectDoppler(false);
  hp->SetProduceFissionFragments(true);
  hp->SetUseWendtFissionModel(false);
  hp->SetUseNRESP71Model(false);

  // Visu
  auto* visManager = new G4VisExecutive();
  visManager->Initialize();

  // UI / Batch
  G4UIExecutive* ui = (argc == 1) ? new G4UIExecutive(argc, argv) : nullptr;
  auto* UImanager = G4UImanager::GetUIpointer();

  auto wall_start = std::chrono::steady_clock::now();
  std::clock_t cpu_start = std::clock();

  if (ui) {
    // Choisis ici le macro par défaut
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
  } else {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);
  }

  auto wall_end = std::chrono::steady_clock::now();
  std::clock_t cpu_end = std::clock();

  const double wall_s = std::chrono::duration<double>(wall_end - wall_start).count();
  const double cpu_s  = double(cpu_end - cpu_start) / CLOCKS_PER_SEC;

  G4cout << "\n==== Timing (process) ====\n"
        << "Wall time : " << wall_s << " s\n"
        << "CPU time  : " << cpu_s  << " s\n"
        << "==========================\n" << G4endl;
  // ---------- SHUTDOWN PROPRE ----------
  // Si tu fermes l’analyse dans les Run/Action via G4AutoDelete, ne fais rien ici.
  // Si, au contraire, tu ouvres/écris/fermes l’analyse dans main, fais-le AVANT de détruire le runManager :
  // {
  //   auto* ana = G4AnalysisManager::Instance();
  //   ana->Write();
  //   ana->CloseFile();
  //   // soit delete ana; soit : G4AutoDelete::Register(ana); (et ne pas le delete ici)
  // }

  delete ui;          // 1) ferme l’UI d’abord
  delete visManager;  // 2) puis la visu
  delete runManager;  // 3) et enfin le run manager (dernier)

  return 0;
}