#pragma once
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
// Needed for using G4ThreeVector in method signatures
#include "G4ThreeVector.hh"

class G4Event;
class G4ParticleGun;
class G4GenericMessenger;

class MyPrimaryGenerator : public G4VUserPrimaryGeneratorAction {
public:
    MyPrimaryGenerator();
    ~MyPrimaryGenerator() override;

    void GeneratePrimaries(G4Event* evt) override;

    // ---------- UI setters (/src/...) ----------
    void SetMode(const G4String& m);            // "mono" | "auto" | "fixed"
    void SetFixedMultiplicity(G4int m) { fFixedMult = (m < 1 ? 1 : m); }
    void SetEqualEnergy(G4bool b)       { fEqualEnergy     = b; }
    void SetUserEnergy(G4double E)      { fUserEnergy_MeV  = E; }        // MeV
    void SetBudget(G4double b)          { fBudget_MeV      = (b > 0 ? b : 30.0); } // MeV
    void SetEnforceBudgetAuto(G4bool b) { fEnforceBudgetAuto = b; }

    // réglages multiplicité (gaussienne)
    void SetNuMean(G4double m)  { fNuMean  = (m > 0 ? m : 1.0); }
    void SetNuSigma(G4double s) { fNuSigma = (s > 0 ? s : 0.5); }

private:
    // gun
    G4ParticleGun*      fParticleGun   = nullptr;
    G4GenericMessenger* fMessenger     = nullptr;

    // Modes
    enum class Mode { kAuto, kMono, kFixed };
    Mode      fMode = Mode::kAuto;

    // Parameters
    G4int     fFixedMult         = 1;      // used in fixed
    G4bool    fEqualEnergy       = false;  // fixed: same E for all?
    G4double  fUserEnergy_MeV    = 2.0;    // mono or fixed+equal
    G4double  fBudget_MeV        = 30.0;   // auto: optional reject if sum>budget
    G4bool    fEnforceBudgetAuto = false;

    // gaussienne pour multiplicité en mode auto
    G4double  fNuMean  = 3.76;   // moyenne Cf-252
    G4double  fNuSigma = 1.2;    // écart-type (ajuste si tu veux)

    // helpers
    G4double       SampleWattMeV();
    G4double SampleWattMeV_Truncated(G4double Emax);
    static G4ThreeVector SampleIsotropicDir();
};