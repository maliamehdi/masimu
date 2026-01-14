#ifndef PhysicsList_h
#define PhysicsList_h

#include "G4VModularPhysicsList.hh"


class MyPhysicsList : public G4VModularPhysicsList
{
public:
	MyPhysicsList();
	~MyPhysicsList() override = default;

public:
	void ConstructParticle() override;
};

#endif
