#include "PrimaryGenerator.hh"
#include "G4GeneralParticleSource.hh"
#include "G4Event.hh"

MyPrimaryGenerator::MyPrimaryGenerator()
{
  fGPS = new G4GeneralParticleSource();
}

MyPrimaryGenerator::~MyPrimaryGenerator()
{
  delete fGPS;
}

void MyPrimaryGenerator::GeneratePrimaries(G4Event* anEvent)
{
  fGPS->GeneratePrimaryVertex(anEvent);
}