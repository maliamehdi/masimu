#pragma once

#include "G4VSensitiveDetector.hh"
#include "G4THitsCollection.hh"
#include "G4VHit.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4SystemOfUnits.hh"

#include <vector>

// ============================
// Hit = 1 cristal (1 copyNo)
// ============================
class CrystalHit : public G4VHit {
public:
  CrystalHit() = default;
  explicit CrystalHit(G4int copy) : fCopyNo(copy) {}

  inline void* operator new(size_t);
  inline void  operator delete(void*);

  void Add(G4double edep, G4double time);

  G4int    GetCopyNo() const { return fCopyNo; }
  G4double GetEdep()   const { return fEdep; }     // MeV
  G4double GetTFirst() const { return fTFirst; }   // ns (ou -1 si aucun dépôt)
  G4double GetTEw()    const { return fTEw; }      // ns (ou -1 si aucun dépôt)

private:
  G4int    fCopyNo  = -1;
  G4double fEdep    = 0.0;   // MeV
  G4double fTFirst  = -1.0;  // ns
  G4double fTEwSum  = 0.0;   // MeV*ns
  G4double fTEw     = -1.0;  // ns
};

using CrystalHitsCollection = G4THitsCollection<CrystalHit>;

extern G4ThreadLocal G4Allocator<CrystalHit>* CrystalHitAllocator;

// ============================
// Sensitive Detector
// ============================
class CrystalSD : public G4VSensitiveDetector {
public:
  CrystalSD(const G4String& name,
            const G4String& hitsCollectionName,
            G4int copyDepth = 0); // depth=0 -> copyNo du volume courant

  ~CrystalSD() override = default;

  void Initialize(G4HCofThisEvent* hce) override;
  G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
  void EndOfEvent(G4HCofThisEvent* /*hce*/) override;

private:
  CrystalHit* GetOrCreateHit(G4int copyNo);

  CrystalHitsCollection* fHitsCollection = nullptr;
  G4int fHCID = -1;
  G4int fCopyDepth = 0;
};