#include "PrimaryGenerator.hh"

MyPrimaryGenerator::MyPrimaryGenerator()
{   
    fParticleGun  = new G4ParticleGun(1);
    //fParticleGun = new G4ParticleSource();

    
    //Neutron source
    
    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition *particle = particleTable->FindParticle("geantino");
    G4ParticleDefinition *gammaray = particleTable->FindParticle("gamma");
    fParticleGun->SetParticleDefinition(particle);
    //fParticleGun->SetParticleDefinition(gammaray);
    
    fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));    
    fParticleGun->SetParticleEnergy(0.*MeV); // initially 0 MeV que j'ai changé
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 0.));
}

MyPrimaryGenerator::~MyPrimaryGenerator()
{
	delete fParticleGun;
}

void MyPrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{
    G4double cosTheta = 2*G4UniformRand() - 1.;
    G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
    G4double phi = CLHEP::twopi*G4UniformRand();
    G4double vx = sinTheta*std::cos(phi),
    vy = sinTheta*std::sin(phi),
    vz = cosTheta;
    
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(vx,vy,vz));
    	
    //Cf 252 source

    if(fParticleGun->GetParticleDefinition() == G4Geantino::Geantino()) 
    {
        G4int Z = 98, A = 252;
        G4double ionCharge   = 0.*eplus;
        G4double excitEnergy = 0.*keV;

        G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
        fParticleGun->SetParticleDefinition(ion);
        fParticleGun->SetParticleCharge(ionCharge);
    }
    else
    {
        fParticleGun->SetParticleEnergy(0.662*MeV); // Cs137
        //fParticleGun->SetParticleEnergy(1.1732*MeV); // Co60
        //fParticleGun->SetParticleEnergy(1.3325*MeV); // Co60
        // Eu-152 (Z=63, A=152), état fondamental
        auto* eu152 = G4IonTable::GetIonTable()->GetIon(63, 152, 0.*keV);
        fParticleGun->SetParticleDefinition(eu152);
        fParticleGun->SetParticleCharge(0.*eplus);

        // direction arbitraire; au repos, RDM gère l’isotropie des produits
        fParticleGun->SetParticleMomentumDirection({0,0,1});
    }

    fParticleGun->GeneratePrimaryVertex(anEvent);
}
