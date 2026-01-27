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
#include <cstdlib>
#include <sys/stat.h>
#include <sys/types.h>

MyRunAction::MyRunAction(const G4String& macroFileName)
: G4UserRunAction(),
  fMacroName(macroFileName)
{
    auto* man = G4AnalysisManager::Instance();

    man->SetVerboseLevel(0); // reduce analysis chatter for speed
    #ifdef G4MULTITHREADED
    man->SetNtupleMerging(true);
    #endif

    // 0) Événements (comptes et énergies agrégées)
    man->CreateNtuple("Events","per-event");
    man->CreateNtupleIColumn("EventID");
    man->CreateNtupleIColumn("NNeutronsEmitted"); // nombre de neutrons émis par la source donc nombre d'évènement * multiplicité aléatoire tiré dans la poissonienne
    man->CreateNtupleDColumn("nDetected"); // nIn nombre de neutrons détectés
    man->CreateNtupleDColumn("EdepCe_keV");
    man->CreateNtupleDColumn("EdepNaI_keV");
    // Hits par ring
    man->CreateNtupleIColumn("HitsRing1");
    man->CreateNtupleIColumn("HitsRing2");
    man->CreateNtupleIColumn("HitsRing3");
    man->CreateNtupleIColumn("HitsRing4");
    // ====== NOUVEAU ======
    man->CreateNtupleDColumn("MeanThermTime_ns");    // moyenne (E<1 eV)
    //man->CreateNtupleDColumn("MeanDetectTime_ns");   // moyenne (entrée ³He), inutile du coup parce que j'enregistre les temps de création des 3H
    man->CreateNtupleIColumn("NNeutronsEscaped");    // compte
    man->CreateNtupleDColumn("lastNeutronTime_ns");    // temps global du dernier neutron (utile pour simuler le deadtime du détecteur)
    
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

  // 5) Neutron primaries: emitted energy and detected ring per neutron
  man->CreateNtuple("NeutronPrimaries", "per primary neutron: E_emit and ring");
  man->CreateNtupleIColumn("eventID");
  man->CreateNtupleIColumn("trackID");
  man->CreateNtupleDColumn("Eemit_MeV");
  man->CreateNtupleIColumn("ring"); // 0 if not detected
  man->CreateNtupleDColumn("detectTime_ns"); // -1 if not detected
  man->FinishNtuple();    // index 5

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

    // 1) Priorité au TAG (fourni par le script bash)
    if (const char* tag = std::getenv("TAG"); tag && *tag) {
        const char* outdirEnv = std::getenv("OUTDIR");
        G4String outDir = "../../myanalyse";
        if (outdirEnv && *outdirEnv) {
            outDir += "/" + G4String(outdirEnv);
        } else {
            outDir += "/" + G4String(tag);
        }

        // mkdir -p outDir
        auto ensureDir = [](const std::string& path){
            if (path.empty()) return;
            std::string cur;
            for (size_t i=0; i<path.size(); ++i) {
                char c = path[i];
                cur.push_back(c);
                if (c=='/' || i==path.size()-1) {
                    if (!cur.empty() && cur!="/") {
                        mkdir(cur.c_str(), 0755);
                    }
                }
            }
        };
        ensureDir(outDir);

        G4String outFile = outDir + "/output_" + G4String(tag) + ".root";
        G4cout << ">>> Ouverture du fichier ROOT (via TAG): " << outFile << G4endl;
        man->OpenFile(outFile);
        return;
    }

    // 2) Fallback: nommage basé sur le macro + runID
    G4String base = fMacroName;               // ex: "run_0.mac"
    if (base.empty()) base = "interactive.mac";
    base = StripPath(base);                   // "run_0.mac"
    base = StripExtension(base, ".mac");      // "run_0"

    std::stringstream tag2;
    tag2 << "_run" << run->GetRunID();        // _run0, _run1, ...

    G4String outFile = "../../myanalyse/FTFP_" + base + tag2.str() + "GATED.root";
    G4cout << ">>> Ouverture du fichier ROOT (fallback): " << outFile << G4endl;
    man->OpenFile(outFile);
}

void MyRunAction::EndOfRunAction(const G4Run*)
{
    auto* man = G4AnalysisManager::Instance();
    man->Write();
    man->CloseFile();
}