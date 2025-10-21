#ifndef SteppingAction_h
#define SteppingAction_h

#include "G4UserSteppingAction.hh"
#include "G4Step.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"

class G4Step;
class MyEventAction; // fwd decl

class MySteppingAction : public G4UserSteppingAction {
public:
  explicit MySteppingAction(MyEventAction* eventAction);
  ~MySteppingAction() override = default;

  void UserSteppingAction(const G4Step*) override;

private:
  MyEventAction* fEventAction = nullptr; // pas utilisé pour l’accumulation d’énergie
};
#endif
