#pragma once

#include "G4UserSteppingAction.hh"
#include "globals.hh"

#include <unordered_set>
#include <unordered_map>

// FWD decl Geant4
class G4Step;
class G4Track;

// FWD decl app
class MyEventAction;

class MySteppingAction : public G4UserSteppingAction {
public:
  explicit MySteppingAction(MyEventAction* eventAction);
  ~MySteppingAction() override = default;

  void UserSteppingAction(const G4Step* step) override;

private:
  // ---------- Helpers "physiques"
  static G4bool IsNeutron(const G4Track* t);
  static G4bool IsBelow1eV(const G4Track* t);
  static G4bool InHe3Cell(const G4Step* step);   // test via tes scoring volumes He-3
  static G4bool AtWorldBoundary(const G4Step* step);

  // ---------- Reset auto au changement d’évènement
  void EnsureEventSync();

  // ---------- États "one-shot" par neutron (indexés par TrackID)
  // 1) Thermalisation : première fois où E_n < 1 eV
  std::unordered_set<G4int>            fThermRecorded;   // déjà vu <1 eV
  std::unordered_map<G4int,G4double>   fThermTime_ns;    // t_therm en ns

  // 2) Détection 3He(n,p)3H : première fois où une capture produit p/t
  std::unordered_set<G4int>            fDetectRecorded;  // capture déjà vue
  std::unordered_map<G4int,G4double>   fDetectTime_ns;   // t_detect en ns

  // 3) Sortie du world : première fois où le neutron franchit la frontière monde
  std::unordered_set<G4int>            fEscapedRecorded; // déjà compté "escaped"

  // ---------- Gestion par-évt
  G4int            fLastEventID = -1;

  // ---------- Lien vers EventAction (pour pousser MeanThermTime, Nescaped, Rings, etc.)
  MyEventAction*   fEventAction = nullptr;
};