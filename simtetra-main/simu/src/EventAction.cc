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

  // Remplissage du ntuple #0 : (eventID, nIn, eCe_keV, eNaI_keV)
  auto* man = G4AnalysisManager::Instance();
  man->FillNtupleIColumn(0, 0, evt->GetEventID());
  man->FillNtupleDColumn(0, 1, nIn);
  man->FillNtupleDColumn(0, 2, eCe/keV);
  man->FillNtupleDColumn(0, 3, eNaI/keV);
  man->AddNtupleRow(0);
  
}