#include "EventAction.hh"
#include "RunAction.hh"    
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"// optionnel si tu ne déréférences pas fRunAction
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4THitsMap.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"
#include "Randomize.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include <unordered_map>
#include <map>
#include <cmath>


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
// ---------- Paramètres de résolution par PARIS ----------
struct ResParams { double resA; double resPower; };

// Valeurs pour le run de calib 152Eu du 05/08/2024 
// static const std::map<std::string, ResParams> parisRes = { //run du 05/08/2024
//     {"PARIS50",  {1.12145,  -0.441244}},
//     {"PARIS70",  {1.80973,  -0.550685}},
//     {"PARIS90",  {1.94868,  -0.564616}},
//     {"PARIS110", {2.11922,  -0.582147}},
//     {"PARIS130", {0.794233, -0.377311}},
//     {"PARIS235", {1.30727,  -0.477402}},
//     {"PARIS262", {1.76345,  -0.542769}},
//     {"PARIS278", {1.98579,  -0.559095}},
//     {"PARIS305", {1.9886,   -0.574021}}
// };

//Avec les copy number
static const std::map<int, ResParams> parisRes = { //run du 05/08/2024
    {0,  {1.12145,  -0.441244}},
    {1,  {1.80973,  -0.550685}},
    {2,  {1.94868,  -0.564616}},
    {3, {2.11922,  -0.582147}},
    {4, {0.794233, -0.377311}},
    {5, {1.30727,  -0.477402}},
    {6, {1.76345,  -0.542769}},
    {7, {1.98579,  -0.559095}},
    {8, {1.9886,   -0.574021}}
};
//Pour les index des différentes partie des PARIS
// ----- Assembly stride/offsets (5 sous-volumes : housing, Ce, NaI, quartz, seal)
static constexpr int kStride   = 1;
static constexpr int kCeOffset = 3;  // Ce = 3,8,13,...
static constexpr int kNaIOffset= 4;  // NaI= 4,7,...

// Renvoie l’index PARIS (0..N-1) à partir du copy number d’un sous-volume donné
inline int ParisIndexFromCopy(int copy, int offset) {
  if (copy < 0) return -1;
  const int d = copy - offset;
  //if (d < 0 || (d % kStride) != 0) return -1;
  return d; // d/kStride;
}
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

// Récupère l'énergie du premier gamma primaire de l'évènement (en keV)
static G4double GetPrimaryGammaEnergyKeV(const G4Event* evt) {
  if (!evt) return 0.0;
  for (G4int iv = 0; iv < evt->GetNumberOfPrimaryVertex(); ++iv) {
    auto* vtx = evt->GetPrimaryVertex(iv);
    if (!vtx) continue;
    for (G4int ip = 0; ip < vtx->GetNumberOfParticle(); ++ip) {
      auto* pp = vtx->GetPrimary(ip);
      if (!pp) continue;
      auto* def = pp->GetParticleDefinition();
      if (def && def->GetParticleName() == "gamma") {
        return pp->GetKineticEnergy()/keV; // → keV
      }
    }
  }
  return 0.0;
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

   // ===== Par-PARIS (imprint) : lire les hits maps & écrire dans le ntuple #3 =====
  auto* hcevt = evt->GetHCofThisEvent();
  //if (!hcevt) goto FILL_SUM_ONLY;

  auto* hmCe  = static_cast<G4THitsMap<G4double>*>(hcevt->GetHC(fHCID_CeEdep));
  auto* hmNaI = static_cast<G4THitsMap<G4double>*>(hcevt->GetHC(fHCID_NaIEdep));

  // Accumulateurs : parisIndex -> (eCe, eNaI) en keV (on convertit directement)
  std::unordered_map<int, std::pair<G4double,G4double>> byParisIndex;

  if (hmCe) {
    for (const auto& kv : *hmCe->GetMap()) {
      const int copy = kv.first;
      //G4cout << "DEBUG: copy number Ce = " << copy-3 << G4endl;
      const int idx  = copy-3;//ParisIndexFromCopy(copy, kCeOffset);
      if (idx < 0) continue; // clef inattendue
      const G4double eMeV = (kv.second ? *(kv.second) : 0.);
      byParisIndex[idx].first += eMeV/keV;  // stocke en keV
    }
  }
  if (hmNaI) {
    for (const auto& kv : *hmNaI->GetMap()) {
      const int copy = kv.first;
      //G4cout << "DEBUG: copy number NaI = " << copy-4 << G4endl;
      const int idx  = copy - 4;//ParisIndexFromCopy(copy, kNaIOffset);
      if (idx < 0) continue;
      const G4double eMeV = (kv.second ? *(kv.second) : 0.);
      byParisIndex[idx].second += eMeV/keV; // keV
    }
  }


  // Récupère les labels « PARIS50 », … depuis la géométrie
  const auto* det =static_cast<const MyDetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  auto* man = G4AnalysisManager::Instance();

  // Pour chaque PARIS, applique la résolution (si >0) et remplis le ntuple #3
  for (const auto& it : byParisIndex) {
    const int idx = it.first;
    const std::string& parisName = det->GetParisLabel(idx);   // "PARIS50", ...

    // Énergies (déjà en keV)
    const double Ece_keV  = it.second.first;
    const double Enai_keV = it.second.second;

    // Paramètres de résolution
    auto prm = parisRes.find(idx);
    if (prm == parisRes.end()) {
      G4Exception("MyEventAction::EndOfEventAction","ParamsNotFound", JustWarning,
                  ("Pas de paramètres de résolution pour " + parisName).c_str());
      // même si pas de params, on peut sauver les valeurs non-smear
      // man->FillNtupleIColumn(3, 0, evt->GetEventID());
      // man->FillNtupleIColumn(3, 1, idx);                // index PARIS
      // man->FillNtupleDColumn(3, 2, Ece_keV);            // sans smear
      // man->FillNtupleDColumn(3, 3, Enai_keV);           // sans smear
      // man->AddNtupleRow(3);
      continue;
    }
    const ResParams& P = prm->second;

    // Smearing uniquement si E>0
    double eResCe_keV  = Ece_keV;
    double eResNaI_keV = Enai_keV;

    if (Ece_keV > 0.0) {
      const double fwhm_Ce  = P.resA * std::pow(Ece_keV,  P.resPower);
      const double sigma_Ce    = (fwhm_Ce / 2.35) * Ece_keV;  // keV
      eResCe_keV  = G4RandGauss::shoot(Ece_keV,  std::max(sigma_Ce,  0.0));
    }
    // if (Enai_keV > 0.0) {
    //   const double fwhm_NaI = P.resA * std::pow(Enai_keV, P.resPower);
    //   const double sigma_NaI   = (fwhm_NaI/ 2.35) * Enai_keV; // keV
    //   eResNaI_keV = G4RandGauss::shoot(Enai_keV, std::max(sigma_NaI, 0.0));
    // }

    // Remplissage ntuple #3 : (eventID, parisIndex, eCe_keV_smear, eNaI_keV_smear)
    man->FillNtupleIColumn(3, 0, evt->GetEventID());
    man->FillNtupleIColumn(3, 1, idx);
    if (eResCe_keV>0)man->FillNtupleDColumn(3, 2, eResCe_keV);
    if (eResNaI_keV>0) man->FillNtupleDColumn(3, 3, eResNaI_keV);
    man->AddNtupleRow(3);
    // ======= NEW : Ntuple "resp" (Etrue/Emeas) =======
    // Pour un seul PARIS à chaque fois
    const G4int ntid = fRunAction->TruthRespNtupleId(); // ← pas de chiffre en dur
    if (ntid >= 0) {
      man->FillNtupleIColumn(ntid, 0, evt->GetEventID());
      man->FillNtupleIColumn(ntid, 1, idx);
      man->FillNtupleDColumn(ntid, 2, Etrue_keV_evt);   // Etrue (mono-énergie GPS)
      man->FillNtupleDColumn(ntid, 3, eResCe_keV);      // Emeas (Ce smearing)
      man->FillNtupleDColumn(ntid, 4, Ece_keV);         // dépôt Ce brut
      man->FillNtupleDColumn(ntid, 5, Enai_keV);        // dépôt NaI brut
      man->AddNtupleRow(ntid);
    }
  }
  //FILL_SUM_ONLY: ; // étiquette pour sauter directement au remplissage des totaux si pas de HCEvt
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

