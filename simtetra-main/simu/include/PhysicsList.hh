#ifndef PhysicsList_h
#define PhysicsList_h

#include "G4VModularPhysicsList.hh"

#include "G4EmPenelopePhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

#include "G4IonElasticPhysics.hh"
#include "G4IonPhysicsXS.hh"
#include "GammaNuclearPhysics.hh"
#include "globals.hh"
#include "NeutronHPphysics.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
//#include "BiasedRDPhysics.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"

class MyPhysicsList : public G4VModularPhysicsList
{
public:
	MyPhysicsList();
	~MyPhysicsList() override = default;

public:
	void ConstructParticle() override;
	//void SetCuts() override;
};

#endif
