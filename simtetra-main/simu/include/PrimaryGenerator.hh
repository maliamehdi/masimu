#ifndef PrimaryGenerator_h
#define PrimaryGenerator_h

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"

#include "G4Geantino.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "globals.hh"

class G4GeneralParticleSource;
class G4Event;

class MyPrimaryGenerator : public G4VUserPrimaryGeneratorAction
{
public:
  MyPrimaryGenerator();
  ~MyPrimaryGenerator() override;

  void GeneratePrimaries(G4Event*) override;

private:
  G4GeneralParticleSource* fGPS = nullptr;
};

#endif