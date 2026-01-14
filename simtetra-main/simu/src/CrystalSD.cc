#include "CrystalSD.hh"

#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TouchableHandle.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ios.hh"

#include <regex>
#include <string>
#include <algorithm> // std::min/max

G4ThreadLocal G4Allocator<CrystalHit>* CrystalHitAllocator = nullptr;

inline void* CrystalHit::operator new(size_t) {
  if (!CrystalHitAllocator) CrystalHitAllocator = new G4Allocator<CrystalHit>;
  return (void*)CrystalHitAllocator->MallocSingle();
}

inline void CrystalHit::operator delete(void* hit) {
  CrystalHitAllocator->FreeSingle((CrystalHit*)hit);
}

void CrystalHit::Add(G4double edep, G4double time) {
  // edep en MeV ; time en ns
  if (edep <= 0.) return;

  // première interaction
  if (fTFirst < 0.0 || time < fTFirst) fTFirst = time;

  fEdep   += edep;
  fTEwSum += edep * time;

  if (fEdep > 0.) fTEw = fTEwSum / fEdep;
}

// ------------------------------------------------------------
//  Helper robuste : extraire l’imprint ID (impr_XX) depuis les noms
//  Exemple typique Geant4 assembly imprint :
//    "av_1_impr_12_SCIONIXPWLVCe_pv_0"
//  -> imprintId = 12
// ------------------------------------------------------------
static G4int ExtractImprintId(const G4TouchableHandle& touch)
{
  if (!touch) return -1;

  // On cherche "impr_<digits>" dans les noms des volumes le long de l'historique
  // On scanne plusieurs profondeurs (0..depth-1) car ça peut être à différents niveaux
  static const std::regex re(R"(impr_(\d+))");

  const int depth = touch->GetHistoryDepth();
  for (int d = 0; d < depth; ++d) {
    const auto* pv = touch->GetVolume(d);
    if (!pv) continue;

    const std::string name = pv->GetName();
    std::smatch m;
    if (std::regex_search(name, m, re)) {
      try {
        return static_cast<G4int>(std::stoi(m[1].str()));
      } catch (...) {
        return -1;
      }
    }
  }
  return -1;
}

CrystalSD::CrystalSD(const G4String& name,
                     const G4String& hitsCollectionName,
                     G4int copyDepth)
: G4VSensitiveDetector(name), fCopyDepth(copyDepth)
{
  collectionName.insert(hitsCollectionName);
}

void CrystalSD::Initialize(G4HCofThisEvent* hce) {
  // 1) créer la collection de hits de l’évènement
  fHitsCollection = new CrystalHitsCollection(SensitiveDetectorName, collectionName[0]);

  // 2) obtenir l’ID une fois
  if (fHCID < 0) {
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  }

  // 3) enregistrer dans l’évènement
  hce->AddHitsCollection(fHCID, fHitsCollection);
}

CrystalHit* CrystalSD::GetOrCreateHit(G4int copyNo) {
  // Recherche linéaire (nb de PARIS petit => OK)
  for (size_t i = 0; i < fHitsCollection->GetSize(); ++i) {
    auto* h = (*fHitsCollection)[i];
    if (h && h->GetCopyNo() == copyNo) return h;
  }
  auto* newHit = new CrystalHit(copyNo);
  fHitsCollection->insert(newHit);
  return newHit;
}

G4bool CrystalSD::ProcessHits(G4Step* step, G4TouchableHistory* /*history*/) {
  if (!step) return false;

  const G4double edep = step->GetTotalEnergyDeposit(); // MeV
  if (edep <= 0.) return false;

  // Temps d’interaction : PreStep (début de l’étape)
  const G4double t = step->GetPreStepPoint()->GetGlobalTime(); // ns

  auto touch = step->GetPreStepPoint()->GetTouchableHandle();
  if (!touch) return false;

  // ------------------------------------------------------------
  //  1) Méthode robuste : imprintId depuis le nom "impr_XX"
  // ------------------------------------------------------------
  G4int copyNo = ExtractImprintId(touch);

  // ------------------------------------------------------------
  //  2) Fallback si imprint introuvable :
  //     (ton ancien comportement)
  // ------------------------------------------------------------
  if (copyNo < 0) {
    copyNo = touch->GetCopyNumber(fCopyDepth);
  }

  // (option debug si tu veux vérifier rapidement)
  // if (copyNo < 0) {
  //   G4cout << "[CrystalSD] copyNo introuvable. depth=" << fCopyDepth
  //          << " pv0=" << touch->GetVolume(0)->GetName() << G4endl;
  // }

  auto* hit = GetOrCreateHit(copyNo);
  hit->Add(edep, t);

  return true;
}

void CrystalSD::EndOfEvent(G4HCofThisEvent* /*hce*/) {
  // Rien à faire : les hits sont déjà prêts (edep, tFirst, tEw)
}