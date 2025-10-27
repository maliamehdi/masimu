// RunAction.cc
#include "RunAction.hh"
#include "EventAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4SystemOfUnits.hh"
#include "G4Types.hh"
#include "G4String.hh"

#include <algorithm>
#include <ctime>

MyRunAction::MyRunAction(const G4String& macroFileName)
: G4UserRunAction(),
  fMacroName(macroFileName)
{
    auto* man = G4AnalysisManager::Instance();

    man->SetVerboseLevel(1);
    #ifdef G4MULTITHREADED
    man->SetNtupleMerging(true);
    #endif

    // 0) Événements (comptes et énergies agrégées)
    man->CreateNtuple("Events","per-event");
    man->CreateNtupleIColumn("EventID");
    man->CreateNtupleDColumn("nThermalEnter"); // nIn nombre de neutrons détectés
    man->CreateNtupleDColumn("EdepCe_keV");
    man->CreateNtupleDColumn("EdepNaI_keV");
    // Hits par ring
    man->CreateNtupleIColumn("HitsRing1");
    man->CreateNtupleIColumn("HitsRing2");
    man->CreateNtupleIColumn("HitsRing3");
    man->CreateNtupleIColumn("HitsRing4");
    man->FinishNtuple(); // index 0

    // 1) Hits triton (par entrée dans une cellule)
    man->CreateNtuple("TritonHits","first-step-in-cell");
    man->CreateNtupleIColumn("EventID");
    man->CreateNtupleDColumn("x_mm");
    man->CreateNtupleDColumn("y_mm");
    man->CreateNtupleDColumn("z_mm");
    man->CreateNtupleDColumn("time_ns");
    man->FinishNtuple(); // index 1

    // 2) Anneaux (un par hit triton)
    man->CreateNtuple("Rings","ring index");
    man->CreateNtupleIColumn("EventID");
    man->CreateNtupleIColumn("RingN");
    man->FinishNtuple(); // index 2

    // 3) Edep par détecteur (par copie)
    man->CreateNtuple("ParisEdep", "Edep par copie");
    man->CreateNtupleIColumn("eventID");
    man->CreateNtupleIColumn("copy");     // copy number
    man->CreateNtupleDColumn("eCe_keV");
    man->CreateNtupleDColumn("eNaI_keV");
    man->FinishNtuple(); // index 3

    // 4) Ntuple Etrue/Emeas (matrice de réponse / migration)
    fTruthRespNtupleId = man->CreateNtuple("resp", "Etrue/Emeas per event & PARIS");
    man->CreateNtupleIColumn(fTruthRespNtupleId, "eventID");
    man->CreateNtupleIColumn(fTruthRespNtupleId, "parisIndex");   // 0..8 (mapping interne)
    man->CreateNtupleDColumn(fTruthRespNtupleId, "Etrue_keV");    // énergie primaire γ
    man->CreateNtupleDColumn(fTruthRespNtupleId, "Emeas_keV");    // énergie mesurée (smeared Ce)
    man->CreateNtupleDColumn(fTruthRespNtupleId, "EdepCe_keV");   // dépôt Ce (avant smearing)
    man->CreateNtupleDColumn(fTruthRespNtupleId, "EdepNaI_keV");  // dépôt NaI (optionnel)
    man->FinishNtuple();    // index 4
}

MyRunAction::~MyRunAction() {}

static G4String StripPath(const G4String& s) {
  std::string ss = s;
  auto pos = ss.find_last_of("/\\");
  if (pos != std::string::npos) ss = ss.substr(pos+1);
  return ss;
}

static G4String StripExtension(const G4String& s, const G4String& ext = ".mac") {
  std::string ss = s;
  if (ss.size() >= ext.size()) {
    if (ss.compare(ss.size()-ext.size(), ext.size(), ext) == 0) {
      ss = ss.substr(0, ss.size()-ext.size());
    }
  }
  return ss;
}

void MyRunAction::BeginOfRunAction(const G4Run* run)
{
    auto* man = G4AnalysisManager::Instance();

    // Construit le nom du fichier de sortie à partir du nom du macro
    G4String base = fMacroName;
    if (base.empty()) base = "interactive.mac";
    base = StripPath(base);
    base = StripExtension(base, ".mac");

    // Ajoute _run<id> si tu veux distinguer plusieurs runs avec la même macro
    std::stringstream tag;
    tag << "_run" << run->GetRunID();

    // Construit le chemin final (adapte le répertoire si besoin)
    G4String outFile = "../../myanalyse/output_" + base + tag.str() + "_smeared.root";

    G4cout << ">>> Ouverture du fichier ROOT : " << outFile << G4endl;
    man->OpenFile(outFile);
}

void MyRunAction::EndOfRunAction(const G4Run*)
{
    auto* man = G4AnalysisManager::Instance();
    man->Write();
    man->CloseFile();
}