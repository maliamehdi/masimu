#include "PhysicsList.hh"
#include "NeutronHPphysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "SToGS_LowEnergyEMPhysicsList.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmParameters.hh"

#include "G4EmPenelopePhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

#include "G4IonElasticPhysics.hh"
#include "G4IonPhysicsXS.hh"
#include "GammaNuclearPhysics.hh"
#include "globals.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

MyPhysicsList::MyPhysicsList()
{	
    RegisterPhysics(new G4DecayPhysics());
    RegisterPhysics(new G4RadioactiveDecayPhysics());
    //RegisterPhysics(new G4EmStandardPhysics());
    //RegisterPhysics(new G4EmStandardPhysics_option4()); 
    
    //RegisterPhysics(new G4EmLivermorePhysics());  // EM low-E soignée
    // + dans ton main ou ici :
    auto em = G4EmParameters::Instance();
    em->SetFluo(true); em->SetAuger(true); em->SetAugerCascade(true); em->SetPixe(true);
    //SetDefaultCutValue(1*um);
    RegisterPhysics(new G4EmPenelopePhysics());
    //RegisterPhysics( new NeutronHPphysics("neutronHP"));
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
