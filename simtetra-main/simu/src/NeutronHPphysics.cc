//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file NeutronHPphysics.cc
/// \brief Implementation of the NeutronHPphysics class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "NeutronHPphysics.hh"

#include "G4GenericMessenger.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessTable.hh"

// Processes

#include "G4HadronElasticProcess.hh"
#include "G4ParticleHPElasticData.hh"
#include "G4ParticleHPThermalScatteringData.hh"
#include "G4ParticleHPElastic.hh"
#include "G4ParticleHPThermalScattering.hh"

#include "G4HadronInelasticProcess.hh"
#include "G4ParticleHPInelasticData.hh"
#include "G4ParticleHPInelastic.hh"

#include "G4NeutronCaptureProcess.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4ParticleHPCapture.hh"

#include "G4NeutronFissionProcess.hh"
#include "G4ParticleHPFissionData.hh"
#include "G4ParticleHPFission.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NeutronHPphysics::NeutronHPphysics(const G4String& name)
: G4VPhysicsConstructor(name)
{
  // define commands for this class
  DefineCommands();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NeutronHPphysics::~NeutronHPphysics()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NeutronHPphysics::ConstructProcess()
{
  G4ParticleDefinition* neutron = G4Neutron::Neutron();
  G4ProcessManager* pManager = neutron->GetProcessManager();
   
  // delete all neutron processes if already registered
  //
  G4VProcess* process = 0;
  process = pManager->GetProcess("hadElastic");
  if (process) pManager->RemoveProcess(process);
  //
  process = pManager->GetProcess("neutronInelastic");
  if (process) pManager->RemoveProcess(process);
  //
  process = pManager->GetProcess("nCapture");      
  if (process) pManager->RemoveProcess(process);
  //
  process = pManager->GetProcess("nFission");      
  if (process) pManager->RemoveProcess(process);      
         
 // (re) create process: elastic
auto* process1 = new G4HadronElasticProcess();
pManager->AddDiscreteProcess(process1);

// --- Thermal scattering (S(a,b)) en premier ---
if (fThermal) {
  // Données S(a,b)
  process1->AddDataSet(new G4ParticleHPThermalScatteringData());

  auto* model1b = new G4ParticleHPThermalScattering();
  model1b->SetMinEnergy(0.0*eV);
  model1b->SetMaxEnergy(4.0*eV); // domaine TS ~ 4 eV
  process1->RegisterMe(model1b);
}

// --- HP elastic "général" au-dessus du domaine thermique ---
process1->AddDataSet(new G4ParticleHPElasticData());

auto* model1a = new G4ParticleHPElastic();
if (fThermal) {
  // Pour éviter le chevauchement avec TS
  model1a->SetMinEnergy(4.0*eV);
}
model1a->SetMaxEnergy(20.0*MeV);
process1->RegisterMe(model1a);


  // (re) create process: inelastic
  //
  G4HadronInelasticProcess* process2 =
    new G4HadronInelasticProcess( "neutronInelastic", G4Neutron::Definition() );
  pManager->AddDiscreteProcess(process2);   
  //
  // cross section data set
  G4ParticleHPInelasticData* dataSet2 = new G4ParticleHPInelasticData();
  process2->AddDataSet(dataSet2);                               
  //
  // models
  G4ParticleHPInelastic* model2 = new G4ParticleHPInelastic();
  process2->RegisterMe(model2);    

  // (re) create process: nCapture   
  //
  G4NeutronCaptureProcess* process3 = new G4NeutronCaptureProcess();
  pManager->AddDiscreteProcess(process3);    
  //
  // cross section data set
  G4ParticleHPCaptureData* dataSet3 = new G4ParticleHPCaptureData();
  process3->AddDataSet(dataSet3);
  //
  // models
  G4ParticleHPCapture* model3 = new G4ParticleHPCapture();
  process3->RegisterMe(model3);
   
  // (re) create process: nFission   
  //
  G4NeutronFissionProcess* process4 = new G4NeutronFissionProcess();
  pManager->AddDiscreteProcess(process4);
  //
  // cross section data set
  G4ParticleHPFissionData* dataSet4 = new G4ParticleHPFissionData();
  process4->AddDataSet(dataSet4);                               
  //
  // models
  G4ParticleHPFission* model4 = new G4ParticleHPFission();
  process4->RegisterMe(model4);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NeutronHPphysics::DefineCommands()
{
  // Define /testhadr/phys command directory using generic messenger class
  fMessenger = new G4GenericMessenger(this,
                        "/testhadr/phys/",
                        "physics list commands");

  // thermal scattering command
  auto& thermalCmd
    = fMessenger->DeclareProperty("thermalScattering", fThermal);

  thermalCmd.SetGuidance("set thermal scattering model");
  thermalCmd.SetParameterName("thermal", false);
  thermalCmd.SetStates(G4State_PreInit);  
}

//..oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
