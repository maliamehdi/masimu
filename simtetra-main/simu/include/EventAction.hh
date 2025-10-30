#pragma once

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <vector>

class MyRunAction;

class MyEventAction : public G4UserEventAction {
public:
  explicit MyEventAction(MyRunAction* runAction);
  ~MyEventAction() override = default;

  // G4 hooks
  void BeginOfEventAction(const G4Event* evt) override;
  void EndOfEventAction(const G4Event* evt) override;

  // Compteurs « anneaux » (appelés depuis SteppingAction)
  void AddHitToRing(G4int ringNumber);
  void ResetRingCounters();

  // Poussé par SteppingAction (neutrons)
  void RegisterNeutronThermalization(G4int /*trackID*/, double t_ns) {
    fThermTimes_ns.push_back(t_ns);
  }
  void RegisterNeutronDetection(G4int /*trackID*/, double t_ns) {
    // enregistre le temps de détection (une fois par neutron)
    fDetectTimes_ns.push_back(t_ns);
  }
  void RegisterNeutronEscaped(G4int /*trackID*/) { ++fNescaped; }
  inline void ResetDetectedNeutrons() { fNDetectedNeutrons = 0; }
  inline void IncrementDetectedNeutrons() { ++fNDetectedNeutrons; }
  inline G4int  GetDetectedNeutrons() const { return fNDetectedNeutrons; }
  inline void SetNNeutronsEmitted(G4int n) { fNNeutronsEmitted = n; }
  inline G4int GetNNeutronsEmitted() const { return fNNeutronsEmitted; }

  // Enregistre la création d'un triton (position en mm, temps en ns)
  inline void RegisterTritonBirth(double x_mm, double y_mm, double z_mm, double t_ns) {
    fTritonBirths.emplace_back(TritonBirth{x_mm, y_mm, z_mm, t_ns});
  }

  // Enregistrement (par trackID neutron) de l'énergie initiale (MeV) et du ring de détection
  inline void RegisterNeutronInitEnergy(G4int neutronTrackID, double E_MeV) {
    fNeutronInitE_MeV[neutronTrackID] = E_MeV;
  }
  inline void RegisterNeutronDetectedRing(G4int neutronTrackID, G4int ring) {
    // Ne pas écraser si déjà connu
    if (fNeutronRing.find(neutronTrackID) == fNeutronRing.end()) {
      fNeutronRing[neutronTrackID] = ring;
    }
  }


private:
  // utilitaires internes
  G4double GetHitsMapSum(G4int hcID, const G4Event* evt) const;

  // pointeur RunAction (pour récupérer l’ID dynamique de l’ntuple "resp")
  MyRunAction* fRunAction = nullptr;

  // IDs des hits collections (résolus lazily)
  G4int fHCID_CeEdep  = -1;
  G4int fHCID_NaIEdep = -1;
  G4int fHCID_CellIn  = -1;

  // Compteurs anneaux (par évènement)
  G4int fHitsRing1 = 0;
  G4int fHitsRing2 = 0;
  G4int fHitsRing3 = 0;
  G4int fHitsRing4 = 0;

  // Buffers neutrons (par évènement)
  std::vector<double> fThermTimes_ns; // temps où E_n < 1 eV (par neutron détecté)
  std::vector<double> fDetectTimes_ns; // temps de détection (entrée ³He)
  G4int fNescaped = 0;
  G4int fNDetectedNeutrons = 0;
  G4int fNNeutronsEmitted = 0;

  // Buffer tritons (par évènement)
  struct TritonBirth { double x_mm; double y_mm; double z_mm; double t_ns; };
  std::vector<TritonBirth> fTritonBirths;

  // Buffers énergies/rings par neutron primaire (clé = TrackID du neutron)
  std::unordered_map<G4int, double> fNeutronInitE_MeV; // énergie initiale émise
  std::unordered_map<G4int, G4int>  fNeutronRing;      // ring 1..4, 0 si non détecté
};