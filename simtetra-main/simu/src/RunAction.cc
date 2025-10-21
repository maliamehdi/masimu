#include "RunAction.hh"
#include "EventAction.hh"

MyRunAction::MyRunAction()
{
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    
    man->SetVerboseLevel(1);
    #ifdef G4MULTITHREADED
    man->SetNtupleMerging(true);
    #endif
    
    // 0) Événements (comptes et énergies agrégées)
    man->CreateNtuple("Events","per-event");
    man->CreateNtupleIColumn("EventID");
    man->CreateNtupleDColumn("nThermalEnter"); // nIn
    man->CreateNtupleDColumn("EdepCe_keV");
    man->CreateNtupleDColumn("EdepNaI_keV");
    man->FinishNtuple(); // index 0

    // 1) Hits triton (par entrée dans une cellule)
    man->CreateNtuple("TritonHits","first-step-in-cell");
    man->CreateNtupleIColumn("EventID");
    man->CreateNtupleDColumn("x_mm");
    man->CreateNtupleDColumn("y_mm");
    man->CreateNtupleDColumn("z_mm");
    man->CreateNtupleDColumn("time_ns");
    man->FinishNtuple(); // index 1

    // 2) Anneaux (un par hit triton)
    man->CreateNtuple("Rings","ring index");
    man->CreateNtupleIColumn("EventID");
    man->CreateNtupleIColumn("RingN");
    man->FinishNtuple(); // index 2
    
}

MyRunAction::~MyRunAction()
{}

void MyRunAction::BeginOfRunAction(const G4Run* run)
{
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    
    G4int runID = run->GetRunID();
    
    std::stringstream strRunID;
    strRunID << runID;
    
    man->OpenFile("../../myanalyse/output"+strRunID.str()+".root");
    
}

void MyRunAction::EndOfRunAction(const G4Run*)
{
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    man->Write();
    man->CloseFile();
    //delete man;
}
