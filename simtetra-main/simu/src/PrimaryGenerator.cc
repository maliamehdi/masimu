#include "PrimaryGenerator.hh"

#include "G4RunManager.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Neutron.hh"
#include "G4IonTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "CLHEP/Random/Randomize.h"

#include <vector>
#include <algorithm>
#include <cmath>

// =======================
// -- Paramètres physiques
// =======================

// Multiplicité moyenne Cf-252 (ajuste si besoin)
static const double kMeanNu = 3.76;

// Paramètres du spectre de Watt W(E;a,b) ~ exp(-E/a) * sinh(sqrt(bE))
static const double kWattA_MeV     = 1.0;  // a en MeV (≈1.0–1.2)
static const double kWattB_invMeV  = 2.3;  // b en 1/MeV (≈2.0–2.5)

// Domaine table pour tirage inverse-CDF
static const double kEmax_MeV      = 20.0;
static const int    kWattNgrid     = 20000;

// *** BUDGET ÉNERGIE TOTALE (somme des E_n primaires) ***
static const double kNeutronBudget_MeV = 30.0;

// *** Multiplicité max pour éviter des boucles inutiles ***
static const int    kMultMax = 14; // (au-delà c’est très rare)

// =======================
// -- Utilitaires internes
// =======================
namespace {

// Table CDF Watt (E, CDF(E)) — construite une fois
struct WattTable {
  std::vector<double> E;    // MeV
  std::vector<double> C;    // CDF normalisée [0,1]
  double a = 0.0, b = 0.0;
  bool built = false;
};

static WattTable gWatt;

inline double WattPDF_unnorm(double E, double a, double b) {
  if (E <= 0.0) return 0.0;
  const double s = std::sqrt(std::max(0.0, b*E));
  return std::exp(-E/a) * std::sinh(s);
}

void BuildWattCDF(double a, double b, double Emax, int N) {
  gWatt.E.resize(N);
  gWatt.C.resize(N);
  const double dE = Emax / (N - 1);

  std::vector<double> pdf(N);
  for (int i=0;i<N;i++) {
    const double E = i * dE;
    gWatt.E[i] = E;
    pdf[i] = WattPDF_unnorm(E, a, b);
  }

  gWatt.C[0] = 0.0;
  for (int i=1;i<N;i++) {
    gWatt.C[i] = gWatt.C[i-1] + 0.5*(pdf[i-1]+pdf[i])*dE;
  }
  const double norm = (gWatt.C.back() > 0.0) ? gWatt.C.back() : 1.0;
  for (int i=0;i<N;i++) gWatt.C[i] /= norm;

  gWatt.a = a; gWatt.b = b; gWatt.built = true;
}

double SampleWattMeV() {
  if (!gWatt.built) BuildWattCDF(kWattA_MeV, kWattB_invMeV, kEmax_MeV, kWattNgrid);
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

inline G4ThreeVector SampleIsotropicDir() {
  const double cosT = 2.0*G4UniformRand() - 1.0;
  const double sinT = std::sqrt(std::max(0.0, 1.0 - cosT*cosT));
  const double phi  = CLHEP::twopi * G4UniformRand();
  return G4ThreeVector(sinT*std::cos(phi), sinT*std::sin(phi), cosT);
}

} // namespace

// =======================
// -- Classe utilisateur
// =======================

MyPrimaryGenerator::MyPrimaryGenerator()
{
  // 1 tir / appel — on boucle dans GeneratePrimaries
  fParticleGun = new G4ParticleGun(1);

  // Particule: neutron
  fParticleGun->SetParticleDefinition(G4Neutron::NeutronDefinition());

  // Source au centre (à adapter si besoin)
  fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleEnergy(1.0*MeV);

  if (!gWatt.built) {
    BuildWattCDF(kWattA_MeV, kWattB_invMeV, kEmax_MeV, kWattNgrid);
  }
}

MyPrimaryGenerator::~MyPrimaryGenerator()
{
  delete fParticleGun;
}

void MyPrimaryGenerator::GeneratePrimaries(G4Event* evt)
{
  // Boucle de rejet : on ne génère l’évènement que s’il respecte le budget
  while (true) {
    // --- multiplicité (≥1, ≤kMultMax)
    int mult = CLHEP::RandPoisson::shoot(kMeanNu);
    
    if (mult < 1) mult = 1;
    if (mult > kMultMax) mult = kMultMax;

    // --- tire les énergies, calcule la somme
    std::vector<double> Elist; Elist.reserve(mult);
    double sumE_MeV = 0.0;
    for (int i=0;i<mult;++i) {
      const double En = SampleWattMeV(); // MeV
      Elist.push_back(En);
      // garder le Elist dans mes tuples de sortie dans EventAction.cc
      sumE_MeV += En;
    }

    // --- test budget
    if (sumE_MeV <= kNeutronBudget_MeV) {
      // OK : on génère les neutrons
      for (double En : Elist) {
        const G4ThreeVector dir = SampleIsotropicDir();
        fParticleGun->SetParticleDefinition(G4Neutron::NeutronDefinition());
        fParticleGun->SetParticleEnergy(En * MeV);
        fParticleGun->SetParticleMomentumDirection(dir);
        fParticleGun->GeneratePrimaryVertex(evt);
      }
      break; // évènement accepté → on sort
    }
    // sinon → on rejette et on retente un autre évènement (nouveaux tirages)
  }
}