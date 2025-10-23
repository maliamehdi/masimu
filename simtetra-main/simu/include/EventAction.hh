#ifndef EventAction_h
#define EventAction_h

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Event;
class MyRunAction; // fwd decl

class MyEventAction : public G4UserEventAction {
public:
  explicit MyEventAction(MyRunAction* runAction); // <-- AJOUT
  ~MyEventAction() override = default;

  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction  (const G4Event*) override;
  
  // Méthodes pour compter les hits par ring
  void AddHitToRing(G4int ringNumber);
  void ResetRingCounters();

private:
  // pointeur vers le RunAction si tu en as besoin (ntuple ids, etc.)
  MyRunAction* fRunAction = nullptr;

  // IDs des hits collections (résolus une fois)
  G4int fHCID_CeEdep  = -1;   // "CeSD/eDep"
  G4int fHCID_NaIEdep = -1;   // "NaISD/eDep"
  G4int fHCID_CellIn  = -1;   // "CellSD/nNeutronEnter"
  
  // Compteurs de hits par ring pour l'événement en cours
  G4int fHitsRing1 = 0;
  G4int fHitsRing2 = 0;
  G4int fHitsRing3 = 0;
  G4int fHitsRing4 = 0;

  // helper
  G4double GetHitsMapSum(G4int hcID, const G4Event* evt) const;
};

#endif