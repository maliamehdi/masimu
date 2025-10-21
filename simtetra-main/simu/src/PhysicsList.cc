#include "PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

MyPhysicsList::MyPhysicsList()
{	
    RegisterPhysics(new G4DecayPhysics());
    RegisterPhysics(new G4RadioactiveDecayPhysics());
    // Radioactive decay
    //RegisterPhysics(new BiasedRDPhysics());

    RegisterPhysics(new G4EmPenelopePhysics());
    RegisterPhysics( new NeutronHPphysics("neutronHP"));
    //RegisterPhysics( new G4HadronPhysicsFTFP_BERT(verb));
    //RegisterPhysics(new G4StoppingPhysics());
    RegisterPhysics(new G4IonPhysics());
}

void MyPhysicsList::ConstructParticle()
{
    G4BosonConstructor  pBosonConstructor;
    pBosonConstructor.ConstructParticle();

    G4LeptonConstructor pLeptonConstructor;
    pLeptonConstructor.ConstructParticle();

    G4MesonConstructor pMesonConstructor;
    pMesonConstructor.ConstructParticle();

    G4BaryonConstructor pBaryonConstructor;
    pBaryonConstructor.ConstructParticle();

    G4IonConstructor pIonConstructor;
    pIonConstructor.ConstructParticle();

    G4ShortLivedConstructor pShortLivedConstructor;
    pShortLivedConstructor.ConstructParticle(); 
}
