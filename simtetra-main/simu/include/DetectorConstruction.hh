#ifndef DetectorConstruction_h
#define DetectorConstruction_h

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polyhedra.hh"
#include "G4Sphere.hh"
#include "G4ExtrudedSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"
#include "G4MultiUnion.hh"
#include <vector>
#include <array>
#include <algorithm>
#include "G4TwoVector.hh"
#include "G4GDMLParser.hh"
#include "G4AssemblyVolume.hh"
#include "G4Transform3D.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackCounter.hh"
#include "G4SDParticleWithEnergyFilter.hh"
#include "G4LogicalVolumeStore.hh"

class MyDetectorConstruction : public G4VUserDetectorConstruction
{

public :
	MyDetectorConstruction();
	~MyDetectorConstruction() override;

    void ConstructSDandField() override;
	virtual G4VPhysicalVolume *Construct() override;
    

    // Logical volumes composant un phoswich PARIS
    G4LogicalVolume* lvCe  = nullptr;   // CeBr3
    G4LogicalVolume* lvNaI = nullptr;   // NaI
    
    G4LogicalVolume *GetScoringVolumeOne() const { return fScoringVolumeOne; }
    G4LogicalVolume *GetScoringVolumeTwo() const { return fScoringVolumeTwo; }
    G4LogicalVolume *GetScoringVolumeThree() const { return fScoringVolumeThree; }
    G4LogicalVolume *GetScoringVolumeFour() const { return fScoringVolumeFour; }
    
    G4LogicalVolume *GetScoringVolume() const { return fScoringVolume; }
    
    G4LogicalVolume *GetSPVolume() const { return fSPVolume; }
    
private:
    
    G4LogicalVolume *logicCellOne, *logicCellTwo, *logicCellThree, *logicCellFour;
	
    G4LogicalVolume *fScoringVolumeOne, *fScoringVolumeTwo, *fScoringVolumeThree, *fScoringVolumeFour;

    G4LogicalVolume *logicCell;
	
    G4LogicalVolume *fScoringVolume;
    
    G4LogicalVolume *logicSP, *logicIC, *logicSupport, *logicShell, *logicShell2, *logic_polycase, *logic_polycase2;
    
    G4LogicalVolume *fSPVolume;
    
    G4double Pressure;
    // A la MAtthieu
    G4VPhysicalVolume* phys_shell;
    G4VPhysicalVolume* phys_shell2;
    G4VPhysicalVolume* phys_polycase;
    G4VPhysicalVolume** phys_cell_gas;
    G4VPhysicalVolume** phys_cell_poly;
    G4VPhysicalVolume* phys_gas;
    
    G4VSolid *solidWorld, *solidModPreSub, *solidHole, *solidLightGuide, *solidUni, *solidMod, *plasticOut, *plasticIn, *solidPlastic, *plexiOut, *plexiIn, *solidPlexi, *solidCore, *solidCoreVac, *solidCrown, *solidComposite, *solidCutBox, *solidFull;
    G4VSolid *solidModA, *solidModB, *solidModPreSubA, *solidModPreSubB, *solidHoleA, *solidHoleB, *solidAirGap, *solidLeft, *solidRight, *solid2Left, *solid2Right;
    G4LogicalVolume *logicModA, *logicModB, *logicAirGap, *logicPolyCell, *logicshieldA, *logicshieldB;
    G4VPhysicalVolume *physModA, *physModB, *physAirGap, *physshieldA, *physshieldB;
    G4SubtractionSolid *modLeft, *modRight, *mod2Left, *mod2Right, *shieldLeft, *shieldRight;

    G4LogicalVolume *logicWorld, *logicMod, *logicUni, *logicPlastic, *logicCyl, *logicPlexi, *logicVac;
    G4VPhysicalVolume *physWorld, *physMod, *physUni, *physPlastic, *physPlexi, *physSP, *physCell1, *physCell1bis, *physCell2, *physCell2bis, *physCell2bisbis, *physCell3, *physCell3bis, *physCell3bisbis, *physCell3bisbisbis, *physCell4, *physCell4bis, *physCell4bisbis, *physCell4bisbisbis, *physCell4bisbisbisbis, *physCyl1, *physCyl1bis, *physCyl2, *physCyl2bis, *physCyl2bisbis, *physCyl3, *physCyl3bis, *physCyl3bisbis, *physCyl3bisbisbis, *physCyl4, *physCyl4bis, *physCyl4bisbis, *physCyl4bisbisbis, *physCyl4bisbisbisbis, *physVac, *monoCell, *monoCyl;
    G4Tubs *solidCell, *solidCyl, *solidVac, *holeLeft, *holeRight;
    G4Sphere *solidSP;
    
    G4GenericMessenger *fMessenger;
    
    G4Isotope *he3Iso;
    G4Element *he3;
    G4Material *co2, *gasMix, *realVac, *air, *modMat, *bore, *borPol, *mod, *plasticMat, *steel, *plexi, *boron, *poly;

    
    
};

#endif
