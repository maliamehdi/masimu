#include "EventAction.hh"
#include "RunAction.hh"
#include "DetectorConstruction.hh"

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

// ---------- Paramètres de résolution par PARIS ----------
struct ResParams { double resA; double resPower; };

// Avec les copy-number (ton tableau)
static const std::map<int, ResParams> parisRes = {
    {0,  {1.12145,  -0.441244}},
    {1,  {1.80973,  -0.550685}},
    {2,  {1.94868,  -0.564616}},
    {3,  {2.11922,  -0.582147}},
    {4,  {0.794233, -0.377311}},
    {5,  {1.30727,  -0.477402}},
    {6,  {1.76345,  -0.542769}},
    {7,  {1.98579,  -0.559095}},
    {8,  {1.9886,   -0.574021}}
};

// Pour les gamma prompt du 252Cf je prends la pire résolution dans le temps pour chaque PARIS
// static const std::map<int, ResParams> parisRes = { //run du 05/08/2024
//     {0,  {1.12145,  -0.441244}}, // 05/08 pour PARIS50 
//     {1,  {1.80973,  -0.550685}}, // 05/08 pour PARIS70 je peux aussi mettre 24/09
//     {2,  {1.65564,	-0.53887}}, //20/06 pour PARIS90
//     {3, {1.93549	-0.564697}}, //29/08 pour PARIS110
//     {4, {0.836128,	-0.368968}}, //17/06 pour PARIS130
//     {5, {1.30727,  -0.477402}}, //05/08 pour PARIS235
//     {6, {1.76345,  -0.542769}}, //05/08 pour PARIS262
//     {7, {1.76703,	-0.536899}}, //07/10 pour PARIS278
//     {8, {1.30165,	-0.503266}} //20/06 pour PARIS305
// };
// ----- Assembly stride/offsets (5 sous-volumes : housing, Ce, NaI, quartz, seal)
static constexpr int kCeOffset  = 3;  // Ce = 3,8,13,...
static constexpr int kNaIOffset = 4;  // NaI= 4,9,14,...

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

// ====== ctor ======
MyEventAction::MyEventAction(MyRunAction* runAction)
: fRunAction(runAction) {}

// Résolution "lazy" des IDs des hits collections (1re fois)
void MyEventAction::BeginOfEventAction(const G4Event* /*evt*/) {
  // Réinitialiser les compteurs de ring au début de chaque événement
  ResetRingCounters();
  ResetDetectedNeutrons();
  fNNeutronsEmitted = 0; 


  // Réinitialiser les buffers neutrons (mêmes lifetime que l’évènement)
  fThermTimes_ns.clear();
  fDetectTimes_ns.clear();
  fNescaped = 0;
  // Réinitialiser le buffer des tritons
  fTritonBirths.clear();
  // Réinitialiser les buffers par-neutron
  fNeutronInitE_MeV.clear();
  fNeutronRing.clear();
  fNeutronDetectTime_ns.clear();

  if (fHCID_CeEdep < 0) {
    auto* sdm = G4SDManager::GetSDMpointer();

    fHCID_CeEdep  = sdm->GetCollectionID("CeSD/edep");
    fHCID_NaIEdep = sdm->GetCollectionID("NaISD/edep");
    fHCID_CellIn  = sdm->GetCollectionID("CellSD/nThermalEnter"); // compteur (filtré)

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

void MyEventAction::AddHitToRing(G4int ringNumber) {
  switch(ringNumber) {
    case 1: ++fHitsRing1; break;
    case 2: ++fHitsRing2; break;
    case 3: ++fHitsRing3; break;
    case 4: ++fHitsRing4; break;
    default: break;
  }
}

void MyEventAction::ResetRingCounters() {
  fHitsRing1 = fHitsRing2 = fHitsRing3 = fHitsRing4 = 0;
}

void MyEventAction::EndOfEventAction(const G4Event* evt) {
  // Totaux par évènement (dep/compteurs)
  const G4double eCe   = GetHitsMapSum(fHCID_CeEdep,  evt);  // MeV
  const G4double eNaI  = GetHitsMapSum(fHCID_NaIEdep, evt);  // MeV
  //const G4double nIn   = GetHitsMapSum(fHCID_CellIn,  evt);  // compteur (filtré n<10 MeV ici)
  // --- nIn = nb de neutrons détectés via ³He(n,p)³H (tritons)
  const G4double nIn = static_cast<G4double>(GetDetectedNeutrons());

  // Énergie gamma primaire (utile pour resp; 0 si source neutrons)
  const double Etrue_keV_evt = GetPrimaryGammaEnergyKeV(evt);

  // Par-PARIS (imprint) : lire les hits maps
  auto* hcevt = evt->GetHCofThisEvent();
  auto* hmCe  = hcevt ? static_cast<G4THitsMap<G4double>*>(hcevt->GetHC(fHCID_CeEdep))  : nullptr;
  auto* hmNaI = hcevt ? static_cast<G4THitsMap<G4double>*>(hcevt->GetHC(fHCID_NaIEdep)) : nullptr;

  // Accumulateurs : parisIndex -> (eCe, eNaI) en keV
  std::unordered_map<int, std::pair<G4double,G4double>> byParisIndex;

  if (hmCe) {
    for (const auto& kv : *hmCe->GetMap()) {
      const int idx = kv.first - kCeOffset;
      if (idx < 0) continue;
      const G4double eMeV = (kv.second ? *(kv.second) : 0.);
      byParisIndex[idx].first += eMeV/keV;
    }
  }
  if (hmNaI) {
    for (const auto& kv : *hmNaI->GetMap()) {
      const int idx = kv.first - kNaIOffset;
      if (idx < 0) continue;
      const G4double eMeV = (kv.second ? *(kv.second) : 0.);
      byParisIndex[idx].second += eMeV/keV;
    }
  }

  auto* man = G4AnalysisManager::Instance();

  // ====== 1) Ntuple #0: Events (ordre qui colle à RunAction.cc) ======
  // Colonnes:
  // 0:EventID, 1:nThermalEnter, 2:EdepCe_keV, 3:EdepNaI_keV,
  // 4:HitsRing1, 5:HitsRing2, 6:HitsRing3, 7:HitsRing4,
  // 8:MeanThermTime_ns, 9:NNeutronsEscaped
  double meanTherm_ns = -1.0;
  if (!fThermTimes_ns.empty()) {
    double sum = 0.0;
    for (double x : fThermTimes_ns) sum += x;
    meanTherm_ns = sum / fThermTimes_ns.size();
  }

  // --- Compte des neutrons primaires émis sur cet évènement
int nEmitted = 0;
for (G4int iv = 0; iv < evt->GetNumberOfPrimaryVertex(); ++iv) {
  auto* vtx = evt->GetPrimaryVertex(iv);
  if (!vtx) continue;
  for (G4int ip = 0; ip < vtx->GetNumberOfParticle(); ++ip) {
    auto* pp  = vtx->GetPrimary(ip);
    if (!pp) continue;
    auto* def = pp->GetParticleDefinition();
    if (def && def->GetParticleName() == "neutron") ++nEmitted;
  }
}

  man->FillNtupleIColumn(0, 0, evt->GetEventID());
  man->FillNtupleIColumn(0, 1, nEmitted); // nombre de neutrons émis par la source
  man->FillNtupleDColumn(0, 2, nIn);
  man->FillNtupleDColumn(0, 3, eCe/keV);
  man->FillNtupleDColumn(0, 4, eNaI/keV);
  man->FillNtupleIColumn(0, 5, fHitsRing1);
  man->FillNtupleIColumn(0, 6, fHitsRing2);
  man->FillNtupleIColumn(0, 7, fHitsRing3);
  man->FillNtupleIColumn(0, 8, fHitsRing4);
  man->FillNtupleDColumn(0, 9, meanTherm_ns);
  man->FillNtupleIColumn(0, 10, fNescaped); // Comment je le calcule ?
  // Dernier temps de détection des neutrons de l'évènement
  double lastDetect_ns = -1.0;
  if (!fDetectTimes_ns.empty()) {
    lastDetect_ns = *std::max_element(fDetectTimes_ns.begin(), fDetectTimes_ns.end());
  }
  man->FillNtupleDColumn(0, 11, lastDetect_ns);

  man->AddNtupleRow(0);

  // ====== 1bis) Ntuple #1: TritonHits (temps et coordonnées à la création)
  // Colonnes: 0:EventID, 1:x_mm, 2:y_mm, 3:z_mm, 4:time_ns
  for (const auto& th : fTritonBirths) {
    man->FillNtupleIColumn(1, 0, evt->GetEventID());
    man->FillNtupleDColumn(1, 1, th.x_mm);
    man->FillNtupleDColumn(1, 2, th.y_mm);
    man->FillNtupleDColumn(1, 3, th.z_mm);
    man->FillNtupleDColumn(1, 4, th.t_ns);
    man->AddNtupleRow(1);
  }

  // ====== 2) Ntuple #3: ParisEdep (par copie) ======
  for (const auto& it : byParisIndex) {
    const int idx = it.first;                 // « copy » logique
    const double Ece_keV  = it.second.first;
    const double Enai_keV = it.second.second;

    auto prm = parisRes.find(idx);
    if (prm == parisRes.end()) {
      // même si pas de params, on peut sauver les valeurs non-smear
      man->FillNtupleIColumn(3, 0, evt->GetEventID());
      man->FillNtupleIColumn(3, 1, idx);
      man->FillNtupleDColumn(3, 2, Ece_keV);
      man->FillNtupleDColumn(3, 3, Enai_keV);
      man->AddNtupleRow(3);
      continue;
    }
    const ResParams& P = prm->second;

    // Smearing Ce (NaI laissé brut, comme ton code)
    double eResCe_keV  = Ece_keV;
    double eResNaI_keV = Enai_keV;

    if (Ece_keV > 0.0) {
      const double fwhm_Ce = P.resA * std::pow(Ece_keV, P.resPower);
      const double sigma_Ce = (fwhm_Ce / 2.35) * Ece_keV;
      eResCe_keV  = G4RandGauss::shoot(Ece_keV, std::max(sigma_Ce, 0.0));
    }

    man->FillNtupleIColumn(3, 0, evt->GetEventID());
    man->FillNtupleIColumn(3, 1, idx);
    if (eResCe_keV  > 0) man->FillNtupleDColumn(3, 2, eResCe_keV);
    if (eResNaI_keV > 0) man->FillNtupleDColumn(3, 3, eResNaI_keV);
    man->AddNtupleRow(3);

    // ====== 3) Ntuple "resp" (id dynamique) ======
    const G4int ntid = fRunAction->TruthRespNtupleId();
    if (ntid >= 0) {
      man->FillNtupleIColumn(ntid, 0, evt->GetEventID());
      man->FillNtupleIColumn(ntid, 1, idx);
      man->FillNtupleDColumn(ntid, 2, Etrue_keV_evt);  // 0 si source neutrons
      man->FillNtupleDColumn(ntid, 3, eResCe_keV);
      man->FillNtupleDColumn(ntid, 4, Ece_keV);
      man->FillNtupleDColumn(ntid, 5, Enai_keV);
      man->AddNtupleRow(ntid);
    }
  }

  // ====== 4) Ntuple #5: NeutronPrimaries (énergie émise et ring détecté)
  // Colonnes: 0:EventID, 1:trackID, 2:E_emit_MeV, 3:ring, 4:detectTime_ns
  for (const auto& kv : fNeutronInitE_MeV) {
    const G4int tid = kv.first;
    const double E_MeV = kv.second;
    G4int ring = 0;
    auto itR = fNeutronRing.find(tid);
    if (itR != fNeutronRing.end()) ring = itR->second;
    double detT_ns = -1.0;
    auto itT = fNeutronDetectTime_ns.find(tid);
    if (itT != fNeutronDetectTime_ns.end()) detT_ns = itT->second;

    man->FillNtupleIColumn(5, 0, evt->GetEventID());
    man->FillNtupleIColumn(5, 1, tid);
    man->FillNtupleDColumn(5, 2, E_MeV);
    man->FillNtupleIColumn(5, 3, ring);
    man->FillNtupleDColumn(5, 4, detT_ns);
    man->AddNtupleRow(5);
  }
}