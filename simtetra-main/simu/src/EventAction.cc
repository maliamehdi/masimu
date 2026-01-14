#include "EventAction.hh"
#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "CrystalSD.hh"

#include "G4RunManager.hh"
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

#include <unordered_map>
#include <map>
#include <cmath>
#include <limits>

// Méthodes pour compter les hits par ring
void MyEventAction::AddHitToRing(G4int ringNumber) {
  switch (ringNumber) {
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

// Avec les copy number
static const std::map<int, ResParams> parisRes = { // run du 05/08/2024
    {1,  {1.12145,  -0.441244}},
    {2,  {1.80973,  -0.550685}},
    {3,  {1.94868,  -0.564616}},
    {4,  {2.11922,  -0.582147}},
    {5,  {0.794233, -0.377311}},
    {6,  {1.30727,  -0.477402}},
    {7,  {1.76345,  -0.542769}},
    {8,  {1.98579,  -0.559095}},
    {9,  {1.9886,   -0.574021}}
};

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

namespace {
  struct ParisAcc {
    G4double eCe_keV   = 0.0;
    G4double eNaI_keV  = 0.0;
    G4double tFirstCe_ns  = -1.0; // -1 si aucun dépôt
    G4double tFirstNaI_ns = -1.0;
  };

  inline void UpdateTFirst(G4double& tFirst, G4double tCand) {
    if (tCand < 0.0) return;
    if (tFirst < 0.0 || tCand < tFirst) tFirst = tCand;
  }

  // Somme d'une HitsMap<G4double> (pour CellSD)
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
  ResetRingCounters();

  if (fHCID_CeEdep < 0) {
    auto* sdm = G4SDManager::GetSDMpointer();

    // >>> NOMS CrystalSD (doivent matcher ConstructSDandField) <<<
    fHCID_CeEdep  = sdm->GetCollectionID("CeCrystalSD/CeCrystalHits");
    fHCID_NaIEdep = sdm->GetCollectionID("NaICrystalSD/NaICrystalHits");

    // CellSD inchangé
    fHCID_CellIn  = sdm->GetCollectionID("CellSD/nThermalEnter");

    if (fHCID_CeEdep < 0 || fHCID_NaIEdep < 0 || fHCID_CellIn < 0) {
      G4Exception("MyEventAction::BeginOfEventAction","MissingHCID", JustWarning,
                  "Au moins une hits collection introuvable. Vérifie les noms.");
    }
  }
}

void MyEventAction::EndOfEventAction(const G4Event* evt) {
  auto* man = G4AnalysisManager::Instance();
  auto* hcevt = evt->GetHCofThisEvent();

  // 1) CellSD : compteur neutron (HitsMap)
  G4double nIn = 0.0;
  if (hcevt && fHCID_CellIn >= 0) {
    auto* hmCell = static_cast<G4THitsMap<G4double>*>(hcevt->GetHC(fHCID_CellIn));
    nIn = SumHitsMap(hmCell);
  }

  // 2) Récupérer collections Crystal
  auto* hcCe  = (hcevt && fHCID_CeEdep  >= 0) ? static_cast<CrystalHitsCollection*>(hcevt->GetHC(fHCID_CeEdep))  : nullptr;
  auto* hcNaI = (hcevt && fHCID_NaIEdep >= 0) ? static_cast<CrystalHitsCollection*>(hcevt->GetHC(fHCID_NaIEdep)) : nullptr;

  // Accumulateurs par idx (ton mapping inchangé)
  std::unordered_map<int, ParisAcc> byParisIndex;

  // Totaux évènement pour ntuple #0
  G4double eCe_evt_MeV  = 0.0;
  G4double eNaI_evt_MeV = 0.0;

  // ---- Ce ----
  if (hcCe) {
    for (size_t i = 0; i < hcCe->GetSize(); ++i) {
      const auto* hit = (*hcCe)[i];
      if (!hit) continue;

      const int copy = hit->GetCopyNo();

      // >>> LOGIQUE INDICES INCHANGEE <<<
      // int idx = copy - 3;
      // if (idx == 21) idx = 8;
      // if (idx < 0) continue;
      int idx = copy;
      // DEBUG sécurité
      //std::cout<<"DEBUG Ce copyNo="<<copy<<" idx="<<idx<<std::endl;
      if (idx < 0 || idx > 8) continue;   // sécurité
      const G4double eMeV   = hit->GetEdep();   // MeV
      const G4double tFirst = hit->GetTFirst(); // ns

      eCe_evt_MeV += eMeV;

      byParisIndex[idx].eCe_keV += eMeV/keV;
      UpdateTFirst(byParisIndex[idx].tFirstCe_ns, tFirst);
    }
  }

  // ---- NaI ----
  if (hcNaI) {
    for (size_t i = 0; i < hcNaI->GetSize(); ++i) {
      const auto* hit = (*hcNaI)[i];
      if (!hit) continue;

      const int copy = hit->GetCopyNo();

      // >>> LOGIQUE INDICES INCHANGEE <<<
      // int idx = copy - 4;
      // if (idx == 21) idx = 8;
      // if (idx < 0) continue;
      int idx = copy;
      if (idx < 0 || idx > 8) continue;   // sécurité

      const G4double eMeV   = hit->GetEdep();
      const G4double tFirst = hit->GetTFirst();

      eNaI_evt_MeV += eMeV;

      byParisIndex[idx].eNaI_keV += eMeV/keV;
      UpdateTFirst(byParisIndex[idx].tFirstNaI_ns, tFirst);
    }
  }

  // Énergie primaire gamma (keV)
  const double Etrue_keV_evt = GetPrimaryGammaEnergyKeV(evt);

  // Labels (si tu en as besoin)
  const auto* det = static_cast<const MyDetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  // 3) Remplissage par idx (ntuple #3, #4, #5)
  for (const auto& it : byParisIndex) {
    const int idx = it.first;
    const ParisAcc& A = it.second;

    const double Ece_keV  = A.eCe_keV;
    const double Enai_keV = A.eNaI_keV;

    // ===== Smearing Ce (inchangé) =====
    double eResCe_keV  = Ece_keV;
    double eResNaI_keV = Enai_keV;

    auto prm = parisRes.find(idx);
    if (prm != parisRes.end()) {
      const ResParams& P = prm->second;
      if (Ece_keV > 0.0) {
        const double resolution_Ce = P.resA * std::pow(Ece_keV, P.resPower);
        const double sigma_Ce      = (resolution_Ce / 2.35) * Ece_keV;
        eResCe_keV = G4RandGauss::shoot(Ece_keV, sigma_Ce);
      }
    } else {
      // pas bloquant
      const std::string& parisName = det->GetParisLabel(idx);
      G4Exception("MyEventAction::EndOfEventAction","ParamsNotFound", JustWarning,
                  ("Pas de paramètres de résolution pour " + parisName).c_str());
    }

    // ===== Ntuple #3 : ParisEdep (eventID, copy(idx), eCe_keV, eNaI_keV) =====
    // >>> IMPORTANT : on remplit TOUJOURS (même 0) <<<
    man->FillNtupleIColumn(3, 0, evt->GetEventID());
    man->FillNtupleIColumn(3, 1, idx);
    man->FillNtupleDColumn(3, 2, Ece_keV);
    man->FillNtupleDColumn(3, 3, Enai_keV);
    man->AddNtupleRow(3);

    // ===== Ntuple #4 : resp (inchangé) =====
    const G4int ntResp = fRunAction->TruthRespNtupleId();
    if (ntResp >= 0) {
      man->FillNtupleIColumn(ntResp, 0, evt->GetEventID());
      man->FillNtupleIColumn(ntResp, 1, idx);
      man->FillNtupleDColumn(ntResp, 2, Etrue_keV_evt);
      man->FillNtupleDColumn(ntResp, 3, eResCe_keV);
      man->FillNtupleDColumn(ntResp, 4, Ece_keV);
      man->FillNtupleDColumn(ntResp, 5, Enai_keV);
      man->AddNtupleRow(ntResp);
    }

    // ===== Ntuple #5 : paris_time (eventID, idx, Ece, Enai, tFirstCe, tFirstNaI) =====
    // >>> tFirst = -1 si pas de dépôt (on écrit quand même) <<<
    man->FillNtupleIColumn(5, 0, evt->GetEventID());
    man->FillNtupleIColumn(5, 1, idx);
    man->FillNtupleDColumn(5, 2, Ece_keV);
    man->FillNtupleDColumn(5, 3, Enai_keV);
    man->FillNtupleDColumn(5, 4, A.tFirstCe_ns);
    man->FillNtupleDColumn(5, 5, A.tFirstNaI_ns);
    man->AddNtupleRow(5);
  }

  // 4) Ntuple #0 : Events (totaux par évènement)
  man->FillNtupleIColumn(0, 0, evt->GetEventID());
  man->FillNtupleDColumn(0, 1, nIn);
  man->FillNtupleDColumn(0, 2, eCe_evt_MeV/keV);
  man->FillNtupleDColumn(0, 3, eNaI_evt_MeV/keV);
  man->FillNtupleIColumn(0, 4, fHitsRing1);
  man->FillNtupleIColumn(0, 5, fHitsRing2);
  man->FillNtupleIColumn(0, 6, fHitsRing3);
  man->FillNtupleIColumn(0, 7, fHitsRing4);
  man->AddNtupleRow(0);
}