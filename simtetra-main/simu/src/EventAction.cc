#include "EventAction.hh"
#include "RunAction.hh"                 // optionnel si tu ne déréférences pas fRunAction
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4THitsMap.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"

// Petit utilitaire pour sommer une HitsMap<G4double>
namespace {
  G4double SumHitsMap(const G4THitsMap<G4double>* hm) {
    if (!hm) return 0.;
    G4double s = 0.;
    for (const auto& kv : *hm->GetMap()) {
      if (kv.second) s += *(kv.second);
    }
    return s;
  }
}

// ====== ctor cohérent avec le .hh ======
MyEventAction::MyEventAction(MyRunAction* runAction)
: fRunAction(runAction) {}

// Résolution "lazy" des IDs des hits collections (1re fois)
void MyEventAction::BeginOfEventAction(const G4Event* /*evt*/) {
  // Réinitialiser les compteurs de ring au début de chaque événement
  ResetRingCounters();
  
  if (fHCID_CeEdep < 0) {
    auto* sdm = G4SDManager::GetSDMpointer();

    // IMPORTANT : ces noms doivent matcher EXACTEMENT ce que tu as enregistré
    // dans ConstructSDandField() : new G4PSEnergyDeposit("edep"), etc.
    fHCID_CeEdep  = sdm->GetCollectionID("CeSD/edep");
    fHCID_NaIEdep = sdm->GetCollectionID("NaISD/edep");
    fHCID_CellIn  = sdm->GetCollectionID("CellSD/nThermalEnter"); // adapte si ton nom diffère

    // (Optionnel) avertir si une collection est introuvable
    if (fHCID_CeEdep < 0 || fHCID_NaIEdep < 0 || fHCID_CellIn < 0) {
      G4Exception("MyEventAction::BeginOfEventAction","MissingHCID", JustWarning,
                  "Au moins une hits collection introuvable. Vérifie les noms.");
    }
  }
}

G4double MyEventAction::GetHitsMapSum(G4int hcID, const G4Event* evt) const {
  if (hcID < 0) return 0.;
  auto* hce = evt->GetHCofThisEvent();
  if (!hce)   return 0.;
  auto* hm = static_cast<G4THitsMap<G4double>*>(hce->GetHC(hcID));
  return SumHitsMap(hm);
}

void MyEventAction::EndOfEventAction(const G4Event* evt) {
  // Totaux par évènement
  const G4double eCe   = GetHitsMapSum(fHCID_CeEdep,  evt);  // MeV
  const G4double eNaI  = GetHitsMapSum(fHCID_NaIEdep, evt);  // MeV
  const G4double nIn   = GetHitsMapSum(fHCID_CellIn,  evt);  // compteur

   // ===== Par-copie : lire les hits maps & écrire dans le ntuple #3 =====
  auto* hparis = evt->GetHCofThisEvent();
  if (!hparis) return;

  auto* hmparisCe  = static_cast<G4THitsMap<G4double>*>(hparis->GetHC(fHCID_CeEdep));
  auto* hmparisNaI = static_cast<G4THitsMap<G4double>*>(hparis->GetHC(fHCID_NaIEdep));

  // Accumulateurs : copyNo -> (eCe, eNaI) en keV
  std::unordered_map<G4int, std::pair<G4double,G4double>> byCopy;

  if (hmparisCe) {
    for (const auto& kv : *hmparisCe->GetMap()) {
      const G4int copy = kv.first;
      const G4double e = kv.second ? *(kv.second) : 0.;
      byCopy[copy].first += e;
    }
  }
  if (hmparisNaI) {
    for (const auto& kv : *hmparisNaI->GetMap()) {
      const G4int copy = kv.first;
      const G4double e = kv.second ? *(kv.second) : 0.;
      byCopy[copy].second += e;
    }
  }

  auto* man = G4AnalysisManager::Instance();

  // Option : récupérer un label lisible depuis la géométrie (sinon fallback "PARIS<copy>")
  // Accès aux labels depuis la géométrie

  //const auto* det = static_cast<const MyDetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  for (const auto& it : byCopy) {
    const G4int copy = it.first;
    const G4double eparisCe  = it.second.first; //MeV
    const G4double eparisNaI = it.second.second; //MeV
    //const std::string& name = det->GetParisLabel(copy);
    // Ntuple #3: (eventID, copy, eCe_keV, eNaI_keV)
    man->FillNtupleIColumn(3, 0, evt->GetEventID());
    man->FillNtupleIColumn(3, 1, copy);
    man->FillNtupleDColumn(3, 2, eparisCe/keV);
    man->FillNtupleDColumn(3, 3, eparisNaI/keV);
    man->AddNtupleRow(3);
  }
  



  // Remplissage du ntuple #0 : (eventID, nIn, eCe_keV, eNaI_keV, HitsRing1, HitsRing2, HitsRing3, HitsRing4)
  
  man->FillNtupleIColumn(0, 0, evt->GetEventID());
  man->FillNtupleDColumn(0, 1, nIn);
  man->FillNtupleDColumn(0, 2, eCe/keV);
  man->FillNtupleDColumn(0, 3, eNaI/keV);
  man->FillNtupleIColumn(0, 4, fHitsRing1);
  man->FillNtupleIColumn(0, 5, fHitsRing2);
  man->FillNtupleIColumn(0, 6, fHitsRing3);
  man->FillNtupleIColumn(0, 7, fHitsRing4);
  man->AddNtupleRow(0);
  

}

// Méthodes pour compter les hits par ring
void MyEventAction::AddHitToRing(G4int ringNumber) {
  switch(ringNumber) {
    case 1: fHitsRing1++; break;
    case 2: fHitsRing2++; break;
    case 3: fHitsRing3++; break;
    case 4: fHitsRing4++; break;
    default: break;
  }
}

void MyEventAction::ResetRingCounters() {
  fHitsRing1 = 0;
  fHitsRing2 = 0;
  fHitsRing3 = 0;
  fHitsRing4 = 0;
}