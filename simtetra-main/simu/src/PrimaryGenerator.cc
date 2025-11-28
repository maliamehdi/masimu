#include "PrimaryGenerator.hh"

#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4Neutron.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"

#include "Randomize.hh"              // G4UniformRand
#include "CLHEP/Random/Randomize.h"  // RandGauss::shoot

#include <vector>
#include <algorithm>
#include <cmath>

// =======================
// -- Paramètres physiques
// =======================

// Multiplicité moyenne Cf-252
static const double kMeanNu  = 3.76;
// Écart-type ~RMS de la multiplicité Cf-252 (voir explications)
static const double kSigmaNu = 1.20;

// Spectre de Watt pour 252Cf: f(E) ~ exp(-E/a) * sinh(sqrt(bE))
static const double kWattA_MeV     = 1.025;
static const double kWattB_invMeV  = 2.926;

// Grille pour inverse-CDF
static const double kEmax_MeV      = 20.0;
static const int    kWattNgrid     = 20000;

// Limiter les cas extrêmes
static const int    kMultMax       = 14;

// =======================
// -- Utilitaires internes
// =======================
namespace {
struct WattTable {
  std::vector<double> E;    // MeV
  std::vector<double> C;    // CDF normalisée [0,1]
  double a=0., b=0.;
  bool built=false;
};
static WattTable gWatt;

inline double WattPDF_unnorm(double E, double a, double b) {
  if (E <= 0.0) return 0.0;
  const double s = std::sqrt(std::max(0.0, b*E));
  return std::exp(-E/a) * std::sinh(s);
}
static void BuildWattCDF(double a, double b, double Emax, int N) {
  gWatt.E.resize(N); gWatt.C.resize(N);
  const double dE = Emax / (N - 1);
  std::vector<double> pdf(N);

  for (int i=0; i<N; ++i) {
    const double E = i * dE;
    gWatt.E[i] = E;
    pdf[i] = WattPDF_unnorm(E, a, b);
  }
  gWatt.C[0] = 0.0;
  for (int i=1; i<N; ++i) gWatt.C[i] = gWatt.C[i-1] + 0.5*(pdf[i-1] + pdf[i]) * dE;

  const double norm = (gWatt.C.back() > 0.0) ? gWatt.C.back() : 1.0;
  for (int i=0; i<N; ++i) gWatt.C[i] /= norm;

  gWatt.a = a; gWatt.b = b; gWatt.built = true;
}
} // namespace

// Tirage inverse-CDF dans la table
G4double MyPrimaryGenerator::SampleWattMeV() {
  if (!gWatt.built) {
    BuildWattCDF(kWattA_MeV, kWattB_invMeV, kEmax_MeV, kWattNgrid);
  }
  const double u = G4UniformRand();
  auto it = std::lower_bound(gWatt.C.begin(), gWatt.C.end(), u);
  if (it == gWatt.C.begin()) return gWatt.E.front();
  if (it == gWatt.C.end())   return gWatt.E.back();
  const int i = int(it - gWatt.C.begin());
  const double C0 = gWatt.C[i-1], C1 = gWatt.C[i];
  const double E0 = gWatt.E[i-1], E1 = gWatt.E[i];
  const double t  = (u - C0) / std::max(1e-16, (C1 - C0));
  return E0 + t*(E1 - E0);
}

// Direction isotrope
G4ThreeVector MyPrimaryGenerator::SampleIsotropicDir() {
  const double cosT = 2.0*G4UniformRand() - 1.0;
  const double sinT = std::sqrt(std::max(0.0, 1.0 - cosT*cosT));
  const double phi  = CLHEP::twopi * G4UniformRand();
  return G4ThreeVector(sinT*std::cos(phi), sinT*std::sin(phi), cosT);
}

// =======================
// -- Classe utilisateur
// =======================
MyPrimaryGenerator::MyPrimaryGenerator()
{
  // Gun: 1 particule par GeneratePrimaryVertex (appelé dans notre boucle)
  fParticleGun = new G4ParticleGun(1);
  fParticleGun->SetParticleDefinition(G4Neutron::NeutronDefinition());
  fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleEnergy(2.0*MeV); // utile pour le mode mono

  if (!gWatt.built) {
    BuildWattCDF(kWattA_MeV, kWattB_invMeV, kEmax_MeV, kWattNgrid);
  }

  // ----- UI messenger (/src/...) -----
  fMessenger = new G4GenericMessenger(this, "/src/", "Primary source control");
  fMessenger->DeclareMethod("mode", &MyPrimaryGenerator::SetMode)
    .SetGuidance("mono | auto | fixed").SetParameterName("mode", false);
  fMessenger->DeclareMethod("mult", &MyPrimaryGenerator::SetFixedMultiplicity)
    .SetGuidance("Fixed multiplicity (mode=fixed)").SetParameterName("m", false).SetDefaultValue("1");
  fMessenger->DeclareMethod("equal", &MyPrimaryGenerator::SetEqualEnergy)
    .SetGuidance("Fixed mode: same energy for all neutrons?").SetParameterName("eq", false).SetDefaultValue("false");
  // Use thread-safe API with explicit default unit (MeV) instead of unit category
  fMessenger->DeclareMethodWithUnit("neutronEnergy", "MeV", &MyPrimaryGenerator::SetUserEnergy)
    .SetGuidance("Mono, or fixed+equal per-neutron energy. Default unit: MeV")
    .SetParameterName("E", false).SetDefaultValue("2");
  fMessenger->DeclareMethod("budgetMEV", &MyPrimaryGenerator::SetBudget)
    .SetGuidance("Auto mode (optional): reject events if sum(E)>budget")
    .SetParameterName("B", false).SetDefaultValue("30");
  fMessenger->DeclareMethod("enforceBudgetAuto", &MyPrimaryGenerator::SetEnforceBudgetAuto)
    .SetGuidance("Auto mode: if true, reject events with sum(E)>budgetMEV")
    .SetParameterName("enf", false).SetDefaultValue("false");
}

MyPrimaryGenerator::~MyPrimaryGenerator() {
  delete fMessenger;
  delete fParticleGun;
}

void MyPrimaryGenerator::SetMode(const G4String& m) {
  if (m=="mono")       fMode = Mode::kMono;
  else if (m=="fixed") fMode = Mode::kFixed;
  else                 fMode = Mode::kAuto; // default/fallback
}

void MyPrimaryGenerator::GeneratePrimaries(G4Event* evt)
{
  // ---------- MODE MONO ----------
  if (fMode == Mode::kMono) {
    fParticleGun->SetParticleEnergy(fUserEnergy_MeV*MeV);
    fParticleGun->SetParticleMomentumDirection(SampleIsotropicDir());
    fParticleGun->GeneratePrimaryVertex(evt);
    return;
  }

  // ---------- CHOIX MULTIPLICITÉ ----------
  int mult = 1;
  if (fMode == Mode::kAuto) {
    // Tirage gaussien tronqué (borné à [1..kMultMax])
    const double draw = CLHEP::RandGauss::shoot(kMeanNu, kSigmaNu);
    mult = (int)std::lround(draw);
    if (mult < 1)        mult = 1;
    if (mult > kMultMax) mult = kMultMax;
  } else { // fixed
    mult = std::max(1, std::min(kMultMax, fFixedMult));
  }

  // ---------- CONSTRUCTION DES ÉNERGIES ----------
  std::vector<double> E_MeV(mult, 0.0);

  if (fMode == Mode::kFixed && fEqualEnergy) {
    // Même énergie pour tous (pas de normalisation automatique)
    const double E = std::max(0.0, (double)fUserEnergy_MeV);
    std::fill(E_MeV.begin(), E_MeV.end(), E);
  } else {
    // auto OU fixed+!equal → tirages Watt (Cf-252)
    auto draw_watt = [&](){
      for (int i=0; i<mult; ++i) E_MeV[i] = SampleWattMeV();
    };

    if (fEnforceBudgetAuto && ((fMode==Mode::kAuto) || (fMode==Mode::kFixed && !fEqualEnergy))) {
      // Rejet si dépasse le budget (appliqué maintenant aussi au mode fixed quand !equal)
      int guard = 0;
      while (true) {
        draw_watt();
        double sum = 0.0; for (double e: E_MeV) sum += e;
        if (sum <= fBudget_MeV) break;
        if (++guard > 128) { // éviter une boucle infinie (un peu plus de tentatives)
          G4cout << "[Budget] Echec: sumE=" << sum << " MeV > " << fBudget_MeV
                 << " MeV après " << guard << " tentatives (mode="
                 << (fMode==Mode::kAuto?"auto":"fixed") << ")" << G4endl;
          break;
        }
      }
    } else {
      draw_watt();
      // E_MeV[0] = 2.0; // SampleWattMeV(); --- IGNORE ---
      // E_MeV[1] = 4.0;
    }
  }

  // ---------- TIRS ----------
  for (int i=0; i<mult; ++i) {
    fParticleGun->SetParticleEnergy(E_MeV[i]*MeV);
    fParticleGun->SetParticleMomentumDirection(SampleIsotropicDir());
    fParticleGun->GeneratePrimaryVertex(evt);
  }
}