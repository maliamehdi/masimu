#ifndef ActionInitialization_h
#define ActionInitialization_h

#include "G4VUserActionInitialization.hh"

#include "PrimaryGenerator.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

class MyActionInitialization : public G4VUserActionInitialization
{
public:
	explicit MyActionInitialization(const G4String& macroFileName = "137Cs.mac");
	~MyActionInitialization();

	
	virtual void Build() const;
	virtual void BuildForMaster() const;
private:
	G4String fMacroName = "137Cs.mac"; // nom du macro pour RunAction

};

#endif
