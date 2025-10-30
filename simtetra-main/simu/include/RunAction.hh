#ifndef RunAction_h
#define RunAction_h

#include "G4UserRunAction.hh"
#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "globals.hh"
#include <sstream>
#include <string>

class MyRunAction : public G4UserRunAction
{
public:
    explicit MyRunAction(const G4String& macroFileName = "standard.mac");
    ~MyRunAction();
    
    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
    // → exposer l'ID de l'ntuple "resp"
    inline G4int TruthRespNtupleId() const { return fTruthRespNtupleId; }

private:
// Utilisé pour nommer le fichier ROOT de sortie
  G4String fMacroName;
  G4int fTruthRespNtupleId = -1; // id de l'ntuple (Etrue/Emeas)
};

#endif
