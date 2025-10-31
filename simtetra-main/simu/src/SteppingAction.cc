#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Neutron.hh"
#include "G4Triton.hh"
#include "G4Proton.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>

MySteppingAction::MySteppingAction(MyEventAction* eventAction)
: fEventAction(eventAction) {}

// --------- Helpers ----------
G4bool MySteppingAction::IsNeutron(const G4Track* t) {
  return t && t->GetParticleDefinition() == G4Neutron::NeutronDefinition();
}
G4bool MySteppingAction::IsBelow1eV(const G4Track* t) {
  return t && t->GetKineticEnergy() < 1.0*eV;
}
G4bool MySteppingAction::InHe3Cell(const G4Step* step) {
  const auto* pre = step->GetPreStepPoint();
  const auto  touch = pre->GetTouchableHandle();
  const auto* lv = touch->GetVolume()->GetLogicalVolume();

  const auto* det =
    static_cast<const MyDetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  auto* v1 = det->GetScoringVolumeOne();
  auto* v2 = det->GetScoringVolumeTwo();
  auto* v3 = det->GetScoringVolumeThree();
  auto* v4 = det->GetScoringVolumeFour();

  return (lv == v1 || lv == v2 || lv == v3 || lv == v4);
}
G4bool MySteppingAction::AtWorldBoundary(const G4Step* step) {
  auto status = step->GetPostStepPoint()->GetStepStatus();
  if (status == fWorldBoundary) return true;
  const auto* postPV = step->GetPostStepPoint()->GetPhysicalVolume();
  return (postPV == nullptr);
}

void MySteppingAction::EnsureEventSync() {
  const auto* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
  const G4int eid = evt ? evt->GetEventID() : -1;
  if (eid != fLastEventID) {
    fLastEventID = eid;
    fThermRecorded.clear();   fThermTime_ns.clear();
    fDetectRecorded.clear();  fDetectTime_ns.clear();
    fEscapedRecorded.clear();
  }
}

// --------- Coeur ----------
void MySteppingAction::UserSteppingAction(const G4Step* step)
{
  EnsureEventSync();

  // --- 1) TEMPS DE THERMALISATION : première fois E_n < 1 eV
  {
    const auto* trk = step->GetTrack();
    if (IsNeutron(trk)) {
      // Enregistrer l'énergie initiale UNIQUEMENT pour les primaires (émis par la source)
      if (trk->GetCurrentStepNumber() == 1 && trk->GetParentID() == 0 && fEventAction) {
        const G4int tid0 = trk->GetTrackID();
        const double Eemit_MeV = trk->GetVertexKineticEnergy()/MeV;
        fEventAction->RegisterNeutronInitEnergy(tid0, Eemit_MeV);
      }

      const G4int tid = trk->GetTrackID();
      if (fThermRecorded.find(tid) == fThermRecorded.end() && IsBelow1eV(trk)) {
        fThermRecorded.insert(tid);
        const double t_ns = step->GetPostStepPoint()->GetGlobalTime()/ns;
        fThermTime_ns[tid] = t_ns;
        if (fEventAction) fEventAction->RegisterNeutronThermalization(tid, t_ns);
      }
    }
  }

  // --- 2) DÉTECTION ³He(n,p)³H : triton/proton créé dans une cellule He-3
  if (InHe3Cell(step)) {
    const auto& secs = *step->GetSecondaryInCurrentStep();
    for (auto* s : secs) {
      auto* def = s->GetDefinition();
      // a) Comptage de détection par neutron parent (une seule fois)
      if (def == G4Triton::TritonDefinition() || def == G4Proton::ProtonDefinition()) {
        const G4int parentNeutronID = s->GetParentID();
        if (parentNeutronID > 0 && fDetectRecorded.find(parentNeutronID) == fDetectRecorded.end()) {
          fDetectRecorded.insert(parentNeutronID);
          const double t_ns = step->GetPostStepPoint()->GetGlobalTime()/ns;
          fDetectTime_ns[parentNeutronID] = t_ns;
          if (fEventAction) {
            fEventAction->RegisterNeutronDetection(parentNeutronID, t_ns);
            fEventAction->IncrementDetectedNeutrons();
          }

        }
      }

      // b) Enregistrement local du triton (la sauvegarde ROOT se fera en fin d'évènement)
      if (def == G4Triton::TritonDefinition()) {
        const auto pos = s->GetPosition();
        const double t_ns = s->GetGlobalTime()/ns;
        if (fEventAction) {
          fEventAction->RegisterTritonBirth(pos.x()/mm, pos.y()/mm, pos.z()/mm, t_ns);
        }
      }
    }
  }

  // --- 3) NEUTRONS QUI SORTENT DU MONDE
  {
    const auto* trk = step->GetTrack();
    if (IsNeutron(trk)) {
      const G4int tid = trk->GetTrackID();
      if (fEscapedRecorded.find(tid) == fEscapedRecorded.end() && AtWorldBoundary(step)) {
        fEscapedRecorded.insert(tid);
        if (fEventAction) fEventAction->RegisterNeutronEscaped(tid);
        // Option: trk->SetTrackStatus(fStopAndKill);
      }
    }
  }

  // --- 4) COMPTAGE DES RINGS (premier pas d’un triton dans une cellule He-3)
  {
    const auto* trk  = step->GetTrack();
    const auto* part = trk->GetParticleDefinition();
    if (part->GetParticleName() == "triton" && step->IsFirstStepInVolume()) {
      const auto* pre = step->GetPreStepPoint();
      auto* lv = pre->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

      const auto* det =
        static_cast<const MyDetectorConstruction*>(
          G4RunManager::GetRunManager()->GetUserDetectorConstruction());

      int ring = 0;
      if      (lv == det->GetScoringVolumeOne())   ring = 1;
      else if (lv == det->GetScoringVolumeTwo())   ring = 2;
      else if (lv == det->GetScoringVolumeThree()) ring = 3;
      else if (lv == det->GetScoringVolumeFour())  ring = 4;

      if (ring > 0 && fEventAction) {
        fEventAction->AddHitToRing(ring);
        // Associer le ring au neutron parent du triton
        const G4int parentNeutronID = trk->GetParentID();
        if (parentNeutronID > 0) {
          fEventAction->RegisterNeutronDetectedRing(parentNeutronID, ring);
        }
      }
    }
  }
}