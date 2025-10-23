#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

MySteppingAction::MySteppingAction(MyEventAction* eventAction)
: fEventAction(eventAction) {}

void MySteppingAction::UserSteppingAction(const G4Step* step)
{
  // Filtrer : seulement si on est dans une des cellules
  const auto* pre = step->GetPreStepPoint();
  auto touch = pre->GetTouchableHandle();
  auto* lv = touch->GetVolume()->GetLogicalVolume();

  const auto* det =
    static_cast<const MyDetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  auto* v1 = det->GetScoringVolumeOne();
  auto* v2 = det->GetScoringVolumeTwo();
  auto* v3 = det->GetScoringVolumeThree();
  auto* v4 = det->GetScoringVolumeFour();
  if (lv != v1 && lv != v2 && lv != v3 && lv != v4) return;

  // Exemple : enregistrer le premier pas d’un triton entrant
  const auto* trk  = step->GetTrack();
  const auto* part = trk->GetParticleDefinition();
  if (part->GetParticleName() != "triton" || !step->IsFirstStepInVolume()) return;

  auto* man = G4AnalysisManager::Instance();

  // --- Ntuple 0: "Hits" -> [EventID(int), fX, fY, fZ, fTime]
  const auto* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
  const G4int eventID = evt ? evt->GetEventID() : -1;
  const auto posW = pre->GetPosition();                 // position MONDE du pas
  const auto t    = pre->GetGlobalTime();

  // ntuple 1 : TritonHits
    man->FillNtupleIColumn(1, 0, eventID);
    man->FillNtupleDColumn(1, 1, posW.x()/mm);
    man->FillNtupleDColumn(1, 2, posW.y()/mm);
    man->FillNtupleDColumn(1, 3, posW.z()/mm);
    man->FillNtupleDColumn(1, 4, t/ns);
    man->AddNtupleRow(1);

    // ntuple 2 : Rings
    // --- Ntuple 2: "Rings" -> [RingN(double)]
    const G4double r = posW.perp()/mm;  // rayon dans le plan transverse (monde)
    // Calcule un entier (1..4) selon r
    G4int ring = 0;
    if      (r >   0. && r <= 100.) ring = 1;
    else if (r > 100. && r <= 150.) ring = 2;
    else if (r > 150. && r <= 200.) ring = 3;
    else if (r > 200. && r <= 250.) ring = 4;
    
    // Remplir le ntuple Rings (comme avant)
    man->FillNtupleIColumn(2, 0, eventID);
    man->FillNtupleIColumn(2, 1, static_cast<G4int>(ring));
    man->AddNtupleRow(2);
    
    // Nouveau : compter le hit pour ce ring dans EventAction
    if (ring > 0) {
      fEventAction->AddHitToRing(ring);
    }
}

  // IMPORTANT :
  //   L’énergie déposée dans Ce/NaI est mesurée par les SDs (G4PSEnergyDeposit)
  //   et consolidée dans MyEventAction via les hits maps.
  // - Si tu veux logger autre chose par pas (par ex. des diagnostics), tu peux
  //   le faire ici, mais garde la physique « principale » dans les SD.
