#include "ActionInitialization.hh"

MyActionInitialization::MyActionInitialization(const G4String& macroFileName)
: G4VUserActionInitialization(),
  fMacroName(macroFileName)
{}
MyActionInitialization::~MyActionInitialization()
{}


void MyActionInitialization::BuildForMaster() const
{
	MyRunAction *runAction = new MyRunAction(fMacroName);
	SetUserAction(runAction);	
}

void MyActionInitialization::Build() const
{
	MyPrimaryGenerator *generator = new MyPrimaryGenerator();
	SetUserAction(generator);
    
    MyRunAction *runAction = new MyRunAction(fMacroName);
    SetUserAction(runAction);
    
    MyEventAction *eventAction = new MyEventAction(runAction);
    SetUserAction(eventAction);
    
    MySteppingAction *steppingAction = new MySteppingAction(eventAction);
    SetUserAction(steppingAction);
}
