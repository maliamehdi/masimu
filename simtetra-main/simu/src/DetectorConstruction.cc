#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Polyhedra.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GDMLParser.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackCounter.hh"
#include "G4SDParticleWithEnergyFilter.hh"

static void PlaceRingCells(
    G4double radius,
    G4double startPhi,
    G4int nSectors,
    G4LogicalVolume* logicCell,
    G4LogicalVolume* logicCyl,
    G4LogicalVolume* logic_polycase,
    G4LogicalVolume* logic_polycase2,
    G4ThreeVector shell_pos,
    G4ThreeVector shell2_pos,
    G4int& copyIndex)
{
    G4double dphi = 360.*deg / nSectors;

    for (G4int i = 0; i < nSectors; ++i)
    {
        G4double phi = i * dphi + startPhi;
        G4double x = radius * std::cos(phi);
        G4double y = radius * std::sin(phi);
        G4ThreeVector globalPos(x, y, 0);
        G4LogicalVolume* logicMod;
        G4ThreeVector localPos;

        if ((x > (+50./6.) && y > 0) || (x > (-50./6.) && y <= 0))
        {
            logicMod = logic_polycase;
            localPos = globalPos - shell2_pos;
        }
        else
        {
            logicMod = logic_polycase2;
            localPos = globalPos - shell_pos;
        }
        // Séparation par X global : si la cellule est à gauche → polycase2, sinon polycase
        
        // if (x < 0) {
        //     logicMod = logic_polycase2;
        //     localPos = globalPos - shell2_pos;
        // } else {
        //     logicMod = logic_polycase;
        //     localPos = globalPos - shell_pos;
        // }
        if (logicMod->GetSolid()->Inside(localPos) != kInside) {
            G4cout << "!!! Cell OUTSIDE volume !!! Pos: " << localPos
                   << " in " << (logicMod == logic_polycase2 ? "polycase2" : "polycase")
                   << G4endl;
        }

        new G4PVPlacement(0, localPos, logicCell, "physCell", logicMod, false, copyIndex++, false);
        new G4PVPlacement(0, localPos, logicCyl,  "physCyl",  logicMod, false, copyIndex++, false);
    }
}

MyDetectorConstruction::MyDetectorConstruction()
{
    fMessenger = new G4GenericMessenger(this, "/detector/", "Detector Construction");
    fMessenger->DeclareProperty("Pressure", Pressure, "Pressure in gas");
    Pressure = 7*6.24151e+08; // MeV/mm3 (non utilisé ici mais conservé)
}

MyDetectorConstruction::~MyDetectorConstruction() {}

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{

    G4Material* ABS = nullptr;
    
    auto nist = G4NistManager::Instance();

    // ===================== Matériaux =====================
    // Air
    auto air     = nist->FindOrBuildMaterial("G4_AIR");
    // Métaux & autres
    steel        = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    auto aluminium = nist->FindOrBuildMaterial("G4_Al");
    auto boron     = nist->FindOrBuildMaterial("G4_B");
    auto methane   = nist->FindOrBuildMaterial("G4_METHANE");
    auto co2       = nist->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
    auto poly_ref  = nist->FindOrBuildMaterial("G4_POLYETHYLENE"); // référence si besoin

    // Éléments standards
    auto elC    = nist->FindOrBuildElement("C");
    auto elN    = nist->FindOrBuildElement("N");
    auto elH    = nist->FindOrBuildElement("H"); // H "normal" (pour ABS)
    // ÉLÉMENT THERMIQUE POUR PE : IMPORTANT
    auto elH_TS = nist->FindOrBuildElement("G4TS_H_of_Polyethylene",/*warning*/false);
    bool hasTS = (elH_TS != nullptr);
    if (hasTS) {
        // FindOrBuild pour être sûr d’avoir l’instance construite
       // elH_TS = nist->FindOrBuildElement("TS_H_of_Polyethylene");
    } else {
        G4cout << "[WARN] TS_H_of_Polyethylene indisponible. "
                "Bascule vers H standard (pas de S(a,b) sur PE dans ce run)."
            << G4endl;
        elH_TS = elH; // fallback
    }


    // ABS (pas d'hydrogène TS ici)
    {
        G4double density = 1.04*g/cm3;
        G4int ncomponents = 3;
        ABS = new G4Material("ABS", density, ncomponents);
        ABS->AddElement(elC, 85*perCent);
        ABS->AddElement(elH,  9*perCent); // H normal
        ABS->AddElement(elN,  6*perCent);
    }

    // He-3
    G4double atomicMass = 3.016*g/mole;
    auto he3Iso = new G4Isotope("he3Iso", 2, 3, atomicMass);
    auto he3    = new G4Element("he3", "he3", 1);
    he3->AddIsotope(he3Iso, 100.*perCent);

    // Gaz cellules He3 + CO2 (inchangé)
    G4double temperatureCell = 293.*kelvin;
    G4double molar_constant  = CLHEP::Avogadro*CLHEP::k_Boltzmann;

    G4double pressureCellOne   = 2*bar;
    G4double densityCellOne    = (atomicMass*pressureCellOne)/(temperatureCell*molar_constant);

    G4double pressureCellTwo   = 2*bar;
    G4double densityCellTwo    = (atomicMass*pressureCellTwo)/(temperatureCell*molar_constant);

    G4double pressureCellThree = 2*bar;
    G4double densityCellThree  = (atomicMass*pressureCellThree)/(temperatureCell*molar_constant);

    G4double pressureCellFour  = 2*bar;
    G4double densityCellFour   = (atomicMass*pressureCellFour)/(temperatureCell*molar_constant);

    auto gasMixOne = new G4Material("he3GasOne", densityCellOne, 2, kStateGas, temperatureCell, pressureCellOne);
    gasMixOne->AddElement(he3, 99.*perCent);
    gasMixOne->AddMaterial(co2, 1.*perCent);

    auto gasMixTwo = new G4Material("he3GasTwo", densityCellTwo, 2, kStateGas, temperatureCell, pressureCellTwo);
    gasMixTwo->AddElement(he3, 99.*perCent);
    gasMixTwo->AddMaterial(co2, 1.*perCent);

    auto gasMixThree = new G4Material("he3GasThree", densityCellThree, 2, kStateGas, temperatureCell, pressureCellThree);
    gasMixThree->AddElement(he3, 99.*perCent);
    gasMixThree->AddMaterial(co2, 1.*perCent);

    auto gasMixFour = new G4Material("he3GasFour", densityCellFour, 2, kStateGas, temperatureCell, pressureCellFour);
    gasMixFour->AddElement(he3, 99.*perCent);
    gasMixFour->AddMaterial(co2, 1.*perCent);

    // ========= Polyethylene avec S(α,β) =========
    // Ratio (CH2) et température fixée à ~293.6 K pour coller aux libs TS
    {
        G4double densityPE = 0.93*g/cm3; // ta valeur initiale
        mod = new G4Material("Polyethylene_TS", densityPE, 2, kStateSolid, 293.6*kelvin);
        mod->AddElement(elC,    1);      // (CH2)n
        mod->AddElement(elH_TS, 2);      // H thermique spécifique PE
        // mod = new G4Material("Polyethylene_ref", densityPE, 2, kStateSolid, 293.6*kelvin);
        // mod->AddMaterial(poly_ref, 1.0); // juste parce que je ne trouve pas TS_H_of_Polyethylene dans ma version de Geant4
    }

    // Polyéthylène boré : partir du PE_TS pour conserver S(α,β) côté H
    borPol = new G4Material("borPol", 1.005*g/cm3, 2, kStateSolid, 293.6*kelvin);
    borPol->AddMaterial(mod,   95.*perCent);
    borPol->AddMaterial(boron, 5.*perCent);

    // ========= Monde =========
    auto solidWorld = new G4Box("solidWorld", 1.*m, 1.*m, 1.*m);
    auto logicWorld = new G4LogicalVolume(solidWorld, air, "logicWorld");
    auto physWorld  = new G4PVPlacement(0, {}, logicWorld, "physWorld", nullptr, false, 0);

    // ===================== Géométrie TETRA =====================
    const G4double R = 36.6*cm;       // non modifié
    const G4double Rshield = 51.6*cm; // non modifié
    const G4double halfHeight = 25.*cm;
    const G4double offset = 148.5*mm;
    const G4double holeRadius = 6.75*cm;

    G4bool checkOverlaps = false;

    G4int numRZ = 4;
    G4double rpos[]  = {0*mm, 516.*mm, 516*mm, 0*mm};
    G4double zpos[]  = {-350*mm,-350.*mm, 350.*mm,350*mm};
    G4double rpos2[] = {0*mm, 366.*mm, 366*mm, 0*mm};
    G4double zpos2[] = {-250.1*mm,-250.1*mm, 250.1*mm, 250.1*mm};

    auto solid_half_shell   = new G4Polyhedra("solid_shell",   270*deg, 180.*deg, 3, numRZ, rpos,  zpos);
    auto solid_half_shell2  = new G4Polyhedra("solid_shell2",   90*deg, 180.*deg, 3, numRZ, rpos,  zpos);
    auto solid_half_polycase  = new G4Polyhedra("solid_half_polycase", 270*deg, 180.*deg, 3, numRZ, rpos2, zpos2);
    auto solid_half_polycase2 = new G4Polyhedra("solid_polycase2",      90*deg, 180.*deg, 3, numRZ, rpos2, zpos2);

    auto box_addon_shell    = new G4Box("shell_addon", 50./3.*mm, 258*mm, 350*mm);
    auto box_addon_polycase = new G4Box("polycase_addon", 50./3.*mm, (183-16.875)*mm, 250*mm);

    auto solid_shell_init1 = new G4UnionSolid("solid_shell_init1", solid_half_shell, box_addon_shell, nullptr, G4ThreeVector(-50./3./2.*mm, -258*mm, 0));
    auto solid_shell_init  = new G4SubtractionSolid("solid_shell_init", solid_shell_init1, box_addon_shell, nullptr, G4ThreeVector(-50./3./2.*mm, 258*mm, 0));

    auto solid_polycase_init1 = new G4UnionSolid("solid_polycase_init1", solid_half_polycase, box_addon_polycase, nullptr, G4ThreeVector(-50./3./2.*mm, (-(183+16.875))*mm, 0));
    auto solid_polycase_init  = new G4SubtractionSolid("solid_polycase_init", solid_polycase_init1, box_addon_polycase, nullptr, G4ThreeVector(-50./3./2.*mm, (183+16.875)*mm, 0));

    auto solid_shell2_init1 = new G4UnionSolid("solid_shell2_init1", solid_half_shell2, box_addon_shell, nullptr, G4ThreeVector(+50./3./2.*mm, 258*mm, 0));
    auto solid_shell2_init  = new G4SubtractionSolid("solid_shell2_init", solid_shell2_init1, box_addon_shell, nullptr, G4ThreeVector(+50./3./2.*mm, -258*mm, 0));

    auto solid_polycase2_init1 = new G4UnionSolid("solid_polycase2_init1", solid_half_polycase2, box_addon_polycase, nullptr, G4ThreeVector(+50./3./2.*mm, (183+16.875)*mm, 0));
    auto solid_polycase2_init  = new G4SubtractionSolid("solid_polycase2_init", solid_polycase2_init1, box_addon_polycase, nullptr, G4ThreeVector(+50./3./2.*mm, (-(183+16.875))*mm, 0));

    auto solid_hole = new G4Tubs("solid_pipe", 0*mm, 67.5*mm, 520.*mm, 0.*deg, 360.0*deg);
    auto solidLightGuide = new G4Box("solidLightGuide", 2.25*cm, 25.*cm, 5.*cm);
    G4ThreeVector yTrans(0., -20.*cm, 0.);
    auto solid_pipe = new G4UnionSolid("solid_pipe", solid_hole, solidLightGuide, nullptr, yTrans);

    auto solid_shell    = new G4SubtractionSolid("solid_shell",    solid_shell_init,    solid_pipe);
    auto solid_polycase = new G4SubtractionSolid("solid_polycase", solid_polycase_init, solid_pipe);
    auto solid_shell2   = new G4SubtractionSolid("solid_shell2",   solid_shell2_init,   solid_pipe);
    auto solid_polycase2= new G4SubtractionSolid("solid_polycase2",solid_polycase2_init,solid_pipe);

    // Volumes logiques : PE_TS sur polycase, boré sur shell
    auto logic_shell     = new G4LogicalVolume(solid_shell,  borPol, "logic_shell");
    auto logic_polycase  = new G4LogicalVolume(solid_polycase,  mod, "logic_polycase");
    auto logic_shell2    = new G4LogicalVolume(solid_shell2, borPol, "logic_shell2");
    auto logic_polycase2 = new G4LogicalVolume(solid_polycase2, mod, "logic_polycase2");

    auto logicLightGuide = new G4LogicalVolume(solidLightGuide, air, "logicLightGuide");

    // ======= Cellules de gaz (inchangé) =======
    auto solidCyl  = new G4Tubs("CylOut", 15.5*mm, 16.*mm, 250.*mm, 0., 360.);
    auto logicCyl  = new G4LogicalVolume(solidCyl, steel, "logicCyl");
    auto solidCell = new G4Tubs("solidCell", 0., 15.5*mm, 250.*mm, 0., 360.);

    auto logicCellOne   = new G4LogicalVolume(solidCell, gasMixOne,   "logicCellOne");
    auto logicCellTwo   = new G4LogicalVolume(solidCell, gasMixTwo,   "logicCellTwo");
    auto logicCellThree = new G4LogicalVolume(solidCell, gasMixThree, "logicCellThree");
    auto logicCellFour  = new G4LogicalVolume(solidCell, gasMixFour,  "logicCellFour");

    fScoringVolumeOne   = logicCellOne;
    fScoringVolumeTwo   = logicCellTwo;
    fScoringVolumeThree = logicCellThree;
    fScoringVolumeFour  = logicCellFour;

    G4int nsect1 = 6;
    G4double dphi1 = 360.*deg;
    G4double ang1 = dphi1/nsect1;
    (void)ang1; // non utilisé mais conservé
    std::vector<G4ThreeVector> hex(nsect1);

    G4double c = 5.*cm;
    G4int copyIndex = 0;
    G4ThreeVector origin = G4ThreeVector(0,0,0);
    

    //PlaceRingCells(..., logicCellOne, logicCyl, logic_polycase, logic_polycase2,
        //       shell_pos, shell2_pos, copyIndex); j'avais origi, origin avant à la place des shell_pos
    PlaceRingCells(2.*c,         -30.*deg, 6, logicCellOne,   logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
    PlaceRingCells(std::sqrt(3.)*c, 0.*deg, 6, logicCellOne,   logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
    PlaceRingCells(3.*c,         -30.*deg, 6, logicCellTwo,   logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
    PlaceRingCells(std::sqrt(7.)*c,  10.893*deg, 6, logicCellTwo, logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
    PlaceRingCells(std::sqrt(7.)*c, -10.893*deg, 6, logicCellTwo, logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
    PlaceRingCells(4.*c,         -30.*deg, 6, logicCellThree, logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
    PlaceRingCells(2.*std::sqrt(3.)*c, 0.*deg, 6, logicCellThree, logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
    PlaceRingCells(std::sqrt(13.)*c,  16.102*deg, 6, logicCellThree, logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
    PlaceRingCells(std::sqrt(13.)*c, -16.102*deg, 6, logicCellThree, logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
    PlaceRingCells(5.*c,         -30.*deg, 6, logicCellFour,  logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
    PlaceRingCells(std::sqrt(19.)*c,  6.587*deg, 6, logicCellFour,  logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
    PlaceRingCells(std::sqrt(21.)*c, 19.107*deg, 6, logicCellFour,  logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
    PlaceRingCells(std::sqrt(21.)*c,-19.107*deg, 6, logicCellFour,  logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
    PlaceRingCells(std::sqrt(19.)*c, -6.587*deg, 6, logicCellFour,  logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);

    
    // Placements
    G4ThreeVector shell_pos  = G4ThreeVector(+offset,0,0);
    G4ThreeVector shell2_pos = G4ThreeVector(-offset,0,0);
    new G4PVPlacement(0, origin,     logic_polycase,  "Polycase",  logic_shell,  false, 0, checkOverlaps);
    new G4PVPlacement(0, origin,     logic_polycase2, "Polycase",  logic_shell2, false, 0, checkOverlaps);
    new G4PVPlacement(0, shell_pos,  logic_shell,     "Shell",     logicWorld,   false, 0, checkOverlaps);
    new G4PVPlacement(0, shell2_pos, logic_shell2,    "Shell2",    logicWorld,   false, 0, checkOverlaps);

    // ====== Châssis alu (inchangé sauf matériaux déjà définis) ======
    
        // Support châssis NORCAN en aluminium de TETRA
        G4Box* norcanlength = new G4Box("norcanlong", 750.*mm, 40*mm, 40*mm);
        G4Box* norcandepth = new G4Box("norcanshort", 40.*mm, 40*mm, 235*mm);
        G4Box* norcanheight = new G4Box("norcanheight", 40.*mm, 375*mm, 40*mm);
        G4Box* norcandepthrenfo = new G4Box("norcanrenfo", 40.*mm, 40*mm, 195*mm);
        G4Box* norcanlengthrenfo = new G4Box("norcanlengthrenfo", 710.*mm, 40*mm, 40*mm);
        //Translations pour les différentes parties du support
        G4ThreeVector ztransback = G4ThreeVector((750.-40)*mm, 0.*mm, -235.*mm); //translation pour le dos
        G4ThreeVector ztransfront = G4ThreeVector((750.-40)*mm, 0.*mm, 235.*mm); //translation pour le devant
        G4ThreeVector ytransright = G4ThreeVector((750.*2.-80.)*mm, 0.*mm, 0.*mm); //translation pour le renforcement droit du support en haut

        G4ThreeVector trans_renfo_bas = G4ThreeVector(0.*mm, -260.*mm, 0.*mm); //translation pour le dos
        G4ThreeVector trans_renfo_bas_droit = G4ThreeVector((750.*2.-80.)*mm, -260.*mm, 0.*mm); //translation pour le devant
        G4ThreeVector trans_renfo_bas_long = G4ThreeVector(0.*mm, -260.*mm, -235.*mm); //translation pour le dos
        G4ThreeVector trans_renfo_bas_long_droit = G4ThreeVector(0.*mm, -260.*mm, 235.*mm); //translation pour le devant

        G4ThreeVector ytransright_bas = G4ThreeVector((750.*2.-80.)*mm, -260.*mm, 0.*mm);
        G4ThreeVector ytransheight1back = G4ThreeVector(0.*mm, -335.*mm, -235.*mm); 
        G4ThreeVector ytransheight1front = G4ThreeVector(0.*mm, -335.*mm, 235.*mm); //translation pour le renforcement haut du support
        G4ThreeVector ytransheight2back = G4ThreeVector((750.*2.-80.)*mm, -335.*mm, -235.*mm); 
        G4ThreeVector ytransheight2front = G4ThreeVector((750.*2.-80.)*mm, -335.*mm,235.*mm); 

        G4UnionSolid* norcan12 = new G4UnionSolid("norcan12",  norcandepth, norcanlength, nullptr, ztransback);
        G4UnionSolid* norcan123 = new G4UnionSolid("norcan123", norcan12, norcanlength, nullptr, ztransfront);
        G4UnionSolid* norcan1234 = new G4UnionSolid("norcan1234", norcan123, norcandepth, nullptr, ytransright);
        G4UnionSolid* norcan12345 = new G4UnionSolid("norcan12345", norcan1234, norcanheight, nullptr, ytransheight1back);
        G4UnionSolid* norcan123456 = new G4UnionSolid("norcan123456", norcan12345, norcanheight, nullptr, ytransheight1front);
        G4UnionSolid* norcan1234567 = new G4UnionSolid("norcan1234567", norcan123456, norcanheight, nullptr, ytransheight2back);
        G4UnionSolid* norcan12345678 = new G4UnionSolid("norcan12345678", norcan1234567, norcanheight, nullptr, ytransheight2front);
        G4UnionSolid* norcanrenfo_1 = new G4UnionSolid("norcanrenfo_1", norcan12345678, norcanlengthrenfo, nullptr, trans_renfo_bas_long);
        G4UnionSolid* norcanrenfo_2 = new G4UnionSolid("norcanrenfo_2", norcanrenfo_1, norcanlengthrenfo, nullptr, trans_renfo_bas_long_droit);
        G4UnionSolid* norcanrenfo_3 = new G4UnionSolid("norcanrenfo_3", norcanrenfo_2, norcandepthrenfo, nullptr, trans_renfo_bas);
        G4UnionSolid* norcanrenfo_4 = new G4UnionSolid("norcanrenfo_4", norcanrenfo_3, norcandepthrenfo,nullptr, trans_renfo_bas_droit);

        //G4UnionSolid* norcanTETRA = ;
        G4LogicalVolume* logicChassis = new G4LogicalVolume(norcan1234, aluminium, "logicChassis");
        G4PVPlacement* physChassis = new G4PVPlacement(0, G4ThreeVector(-710.*mm, -590.*mm, 0.*mm), logicChassis, "physChassis", logicWorld, false, 0, checkOverlaps);
        // //---------------------------------------------------
        // // Support châssis NORCAN en aluminium de la chambre d'ionisation + PARIS
        // G4Box* barre_1 = new G4Box("barre_1", 30.*mm, 915*mm, 30*mm);//grande barre
        // G4Box* barre_2 = new G4Box("barre_2", 30.*mm, 915*mm, 30*mm); //petite barre barre
        // G4Box* barre_3 = new G4Box("barre_3", 30.*mm, 30*mm, 500*mm); // barre horizontale du dessus
        // G4Box* barre_4 = new G4Box("barre_4", 30.*mm, 400*mm, 30*mm); //renforcement verticalextéreieur
        // G4Box* barre_5 = new G4Box("barre_5", 20.*mm, 300*mm, 20*mm); //barre verticale soutenant la chambre d'ionisation
        // G4Box* barre_6 = new G4Box("barre_6", 30.*mm, 30*mm, 65.05*mm); //petite barre  horizontale soutenant le front de la chambre d'ionisation
        // G4Box* barre_7 = new G4Box("barre_7", 30.*mm, 30*mm, 74.7*mm); //petite barre horizontale soutenant le dos de la chambre d'ionisation
        // G4Box* barre_8 = new G4Box("barre_8", 20.*mm, 20*mm, 95.25*mm); //petite barre horizontale de soutien sous la chambre d'ionisation
        // G4Box* barre_9 = new G4Box("barre_9", 30.*mm, 30*mm, 82.5*mm); //petite barre horizontale de soutien exterieur du support 
        
        // G4UnionSolid* support1 = new G4UnionSolid("support1", barre_3, barre_1, nullptr, G4ThreeVector(0., -945, -305.));
        // G4UnionSolid* support2 = new G4UnionSolid("support2", support1, barre_1, nullptr, G4ThreeVector(0., -945, 305.));
        // G4UnionSolid* support3 = new G4UnionSolid("support3", support2, barre_4, nullptr, G4ThreeVector(0., -430., -470.));
        // G4UnionSolid* support4 = new G4UnionSolid("support4", support3, barre_4, nullptr, G4ThreeVector(0., -430., 470.));
        // G4UnionSolid* support5 = new G4UnionSolid("support5", support4, barre_7, nullptr, G4ThreeVector(0., -840., -219.95));
        // G4UnionSolid* support6 = new G4UnionSolid("support6", support5, barre_6, nullptr, G4ThreeVector(0., -840., 200.3));
        // G4UnionSolid* support7 = new G4UnionSolid("support7", support6, barre_5, nullptr, G4ThreeVector(0., -620., -135.6));
        // G4UnionSolid* support8 = new G4UnionSolid("support8", support7, barre_5, nullptr, G4ThreeVector(0., -620., 174.9));
        
        // G4LogicalVolume* logicSupport = new G4LogicalVolume(support8, aluminium, "logicSupport");
        // // Placement du support
        // G4PVPlacement* physSupport = new G4PVPlacement(0, G4ThreeVector(0.*mm, 590.*mm, 0.*mm), logicSupport, "physSupport", logicWorld, false, 0, checkOverlaps);
        // géométries (demi-dimensions)
        auto s_barre_1 = new G4Box("s_barre_1", 30*mm, 800*mm, 30*mm);   // 60x1600
        auto s_barre_2 = new G4Box("s_barre_2", 30*mm, 800*mm, 30*mm);
        auto s_barre_3 = new G4Box("s_barre_3", 30*mm, 30*mm, 500*mm);   // 60x1000 (traverse)
        auto s_barre_4 = new G4Box("s_barre_4", 30*mm, 400*mm, 30*mm);   // 60x800 (renfort)
        auto s_barre_5 = new G4Box("s_barre_5", 20*mm, 300*mm, 20*mm);  // 40x600 (VERTICALE en Y)
        auto s_barre_6 = new G4Box("s_barre_6", 30*mm, 30*mm, 65.05*mm); // 60x130.1
        auto s_barre_7 = new G4Box("s_barre_7", 30*mm, 30*mm, 74.7*mm);  // 60x149.4
        auto s_barre_8 = new G4Box("s_barre_8", 20*mm, 20*mm, 115.25*mm); // 40x190.5 (sous la chambre)
        auto s_barre_9 = new G4Box("s_barre_9", 30*mm, 30*mm, 82.5*mm);  // 60x165 (exterieur du support)
        // logicals (alu)
        auto lv_barre_1 = new G4LogicalVolume(s_barre_1, aluminium, "lv_barre_1");
        auto lv_barre_2 = new G4LogicalVolume(s_barre_2, aluminium, "lv_barre_2");
        auto lv_barre_3 = new G4LogicalVolume(s_barre_3, aluminium, "lv_barre_3");
        auto lv_barre_4 = new G4LogicalVolume(s_barre_4, aluminium, "lv_barre_4");
        auto lv_barre_5 = new G4LogicalVolume(s_barre_5, aluminium, "lv_barre_5");
        auto lv_barre_6 = new G4LogicalVolume(s_barre_6, aluminium, "lv_barre_6");
        auto lv_barre_7 = new G4LogicalVolume(s_barre_7, aluminium, "lv_barre_7");
        auto lv_barre_8 = new G4LogicalVolume(s_barre_8, aluminium, "lv_barre_8");
        auto lv_barre_9 = new G4LogicalVolume(s_barre_9, aluminium, "lv_barre_9");

        // base = ancien offset du mother : (0, +590, 0) + (115-20, 0, 0) rayon du support chambre ABS
        const G4ThreeVector base(115-33.5*mm, 550.*mm, 0.); //115mm rayon deu support ABS chambre - 33.5mm j'ai suivi le schéma

        // placements (add base.y à toutes tes Y locales)
        new G4PVPlacement(nullptr, base + G4ThreeVector(0,     0*mm,   0*mm),
                        lv_barre_3, "pv_barre_3", logicWorld, false, 0, checkOverlaps);

        new G4PVPlacement(nullptr, base + G4ThreeVector(0,  -830*mm, -305*mm),
                        lv_barre_1, "pv_barre_1", logicWorld, false, 0, checkOverlaps);
        new G4PVPlacement(nullptr, base + G4ThreeVector(0,  -830*mm, +305*mm),
                        lv_barre_2, "pv_barre_2", logicWorld, false, 0, checkOverlaps);

        new G4PVPlacement(nullptr, base + G4ThreeVector(0,  -430*mm, -470*mm),
                        lv_barre_4, "pv_barre_4_L", logicWorld, false, 0, checkOverlaps);
        new G4PVPlacement(nullptr, base + G4ThreeVector(0,  -430*mm, +470*mm),
                        lv_barre_4, "pv_barre_4_R", logicWorld, false, 1, checkOverlaps);

        new G4PVPlacement(nullptr, base + G4ThreeVector(0,  -840*mm, -200.95*mm),
                        lv_barre_7, "pv_barre_7", logicWorld, false, 0, checkOverlaps);
        new G4PVPlacement(nullptr, base + G4ThreeVector(0,  -840*mm, +220.3*mm),
                        lv_barre_6, "pv_barre_6", logicWorld, false, 0, checkOverlaps);

        // barre_5 = montants chambre VERTICAUX (axe Y)
        new G4PVPlacement(nullptr, base + G4ThreeVector(0*mm,  -838*mm, -130.6*mm), 
                        lv_barre_5, "pv_barre_5_L", logicWorld, false, 0, checkOverlaps);
        new G4PVPlacement(nullptr, base + G4ThreeVector(0*mm,  -838*mm, +138.9*mm),
                        lv_barre_5, "pv_barre_5_R", logicWorld, false, 1, checkOverlaps);
        //barre de renforcement sous la chambre
        new G4PVPlacement(nullptr, base + G4ThreeVector(0,  -950*mm, 5.*mm),
                        lv_barre_8, "pv_barre_8", logicWorld, false, 0, checkOverlaps);
        //barre de renforcement exterieure
                new G4PVPlacement(nullptr, base + G4ThreeVector(0,  -860*mm, -417.5*mm),
                        lv_barre_9, "pv_barre_9_L", logicWorld, false, 0, checkOverlaps);
        new G4PVPlacement(nullptr, base + G4ThreeVector(0,  -860*mm, +417.5*mm),
                        lv_barre_9, "pv_barre_9_R", logicWorld, false, 1, checkOverlaps); 
        
        
        // Rail courbé PARIS
        G4double R_arc    = 550.*mm;     // rayon moyen de l’arc
        G4double w_arc    = 50.*mm;      // largeur radiale (dans le plan)
        G4double t    = 20.*mm;      // épaisseur (hors plan, le long de Z)
        G4double phi_arc = 20.*deg;     // angle de départ
        G4double dphi_arc = 140.*deg;     // ouverture de l’arc

        G4double rmin_arc = R_arc - w_arc/2.;
        G4double rmax_arc = R_arc + w_arc/2.;

        auto solidArc = new G4Tubs("Arc",
                                rmin_arc, rmax_arc,
                                t/2.,
                                phi_arc, dphi_arc);

        auto logicArc = new G4LogicalVolume(solidArc, aluminium, "logicArc");
        G4RotationMatrix* rotY = new G4RotationMatrix();
        G4RotationMatrix* rotY2 = new G4RotationMatrix();
        rotY->rotateY(90.*deg);  // fait passer l’axe Z → X
        rotY2->rotateY(90.*deg); // fait passer l’axe Z → X pour l'autre arc
        rotY2->rotateZ(180.*deg); // fait passer l’axe Z → -X pour l'autre arc
        G4PVPlacement* physArc1 = new G4PVPlacement(rotY, G4ThreeVector(-40+81.5*mm, 0.*mm, 0.*mm), logicArc, "physArc1", logicWorld, false, 0, false);
        G4PVPlacement* physArc2 = new G4PVPlacement(rotY2, G4ThreeVector(-40+81.5*mm, 0.*mm, 0.*mm), logicArc, "physArc2", logicWorld, false, 1, false);

        //---------------------------------------------------
        // Support en impression 3D plastique  la chambre d'ionisation ABS  (acrylonitrile butadiène styrène)
        //     Notes rapides
	    // •	Ajuste w_arc pour que l’arc épouse bien ton diamètre de détecteur (par ex. w_arc ≈ diamètre_du_tube si tu veux un vrai « berceau »).
	    // •	Si tu veux que l’arc soit derrière au lieu de devant, mets la translation à -dZ.
	    // •	Si l’arc ne doit pas couvrir exactement 180°, change dphi.
        // --- Paramètres ---
        const G4double Rmean_supp    = 115.*mm;  // rayon moyen demandé
        const G4double t_supp        = 24.*mm;   // épaisseur (au total) le long de l’axe
        const G4double phi0_supp     = 180.*deg; // 180° → 360°
        const G4double dphi_supp     = 180.*deg;

        // Largeur radiale de l’arc d’extension (à choisir selon ton détecteur)
        const G4double w_arc_supp    = 30.*mm;   // 

        // Demi-épaisseurs pour G4Tubs
        const G4double hz       = t_supp/2.;     // 12 mm
        const G4double hz_ext   = t_supp/2.;     // extension avec même épaisseur

        // --- Solide du demi-cylindre plein (berceau principal) ---
        auto solidSemi = new G4Tubs("SupportSemi",
                                    /*rmin*/ 0.*mm,
                                    /*rmax*/ Rmean_supp,
                                    /*hz*/   hz,
                                    /*phi*/  phi0_supp, dphi_supp);

        auto logicSemi = new G4LogicalVolume(solidSemi, ABS, "logicSupportABS");

        // --- Arc d’extension : secteur annulaire centré sur Rmean ---
        // rmin / rmax définis pour avoir rayon moyen = Rmean et largeur radiale = w_arc
        const G4double rmin_ext_supp = Rmean_supp - w_arc_supp/2.;
        const G4double rmax_ext_supp = Rmean_supp + w_arc_supp/2.;

        auto solidArcExt = new G4Tubs("SupportArcExt",
                                    rmin_ext_supp, rmax_ext_supp,
                                    hz_ext,
                                    phi0_supp, dphi_supp);

        // --- Union : on colle l’arc EXT devant le demi-cylindre ---
        // Translation axiale = + (hz + hz_ext) pour le placer juste au contact
        G4ThreeVector dZ(0,0, hz + hz_ext);

        auto solidSupportWithArc =
        new G4UnionSolid("SupportABS_WithArc", solidSemi, solidArcExt, nullptr, dZ);
        // ---------- TROU demi-cylindre ----------
        const G4double R_hole      = 70.*mm;  // rayon du trou
        const G4double offsetHoleY = 25.*mm;  // 25 mm "en dessous" (vers -Y)

        // On prend un demi-cylindre (même secteur angulaire que le support) et
        // on le fait un poil plus long que la somme des épaisseurs pour traverser proprement.
        const G4double hz_hole = (t/2.) + (t/2.) + 1.*mm;  // un peu plus que hz + hz_ext

        auto solidHole = new G4Tubs("SupportHoleSemi",
                                    /*rmin*/ 0.*mm,
                                    /*rmax*/ R_hole,
                                    /*hz*/   hz_hole,
                                    /*phi*/  phi0_supp, dphi_supp);  // 180°→360° comme le support

        // Soustraction : translation locale (0, -25 mm, 0) = vers le bas
        auto solidSupportWithArc_Hole =
        new G4SubtractionSolid("SupportABS_WithArc_Hole",
                                solidSupportWithArc,
                                solidHole,
                                /*rot*/ nullptr,
                                /*trans*/ G4ThreeVector(0., -offsetHoleY, 0.));

        // (re)crée le logical avec le solide percé
        auto logicSupportWithArc =
        new G4LogicalVolume(solidSupportWithArc_Hole, ABS, "logicSupportABS_WithArc");
        

        // --- Placements d’exemple (axe = Z par défaut) ---
        new G4PVPlacement(nullptr, G4ThreeVector(0,0,-100*mm),
                        logicSupportWithArc, "physSupportABS_A",
                        logicWorld, false, 0, true);
        // On fait une rotation de 180° pour le second placement selon X
        G4RotationMatrix* flip = new G4RotationMatrix();
        flip->rotateX(180.*deg);  // rotation de 180° autour de l'axe X
        flip->rotateZ(180.*deg); // rotation de 180° autour de l'axe Y pour le second placement
        // Placement du second support
        new G4PVPlacement(flip, G4ThreeVector(0,0,+100*mm),
                        logicSupportWithArc, "physSupportABS_B",
                        logicWorld, false, 1, true);
    // ====== Chambre d’ionisation (inchangé) ======
    G4double rICint = 89.5*mm;
    G4double rICext = 90.*mm;
    G4double hIC    = 77.675*mm;

    auto solidIC    = new G4Tubs("solidIC", rICint, rICext, hIC, 0., 360.*deg);
    auto solidICend = new G4Tubs("solidICend", 0., rICext+10*mm, 13.37*mm, 0., 360.*deg);
    auto solidICfinal = new G4UnionSolid("solidICfinal", solidIC, solidICend, nullptr, G4ThreeVector(0., 0., hIC + 13.37*mm));
    auto logicIC = new G4LogicalVolume(solidICfinal, aluminium, "logicIC");

    auto solidMethane = new G4Tubs("solidMethane", 0., rICint, hIC, 0., 360.*deg);
    auto logicMethane = new G4LogicalVolume(solidMethane, methane, "logicMethane");
    (void)new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicMethane, "physMethane", logicIC, false, 0, checkOverlaps);
    (void)new G4PVPlacement(0, origin, logicIC, "physIC", logicWorld, false, 0, checkOverlaps);

    // ====== Source Point (inchangé) ======
    solidSP  = new G4Sphere("SP", 0., 0.1*mm, 0., 360., 0., 180.);
    logicSP  = new G4LogicalVolume(solidSP, air, "logicSP");
    physSP   = new G4PVPlacement(0, G4ThreeVector(0., -2.*cm, 10.*cm), logicSP, "physSP", logicWorld, false, 0, true);
    fSPVolume = logicSP;

    // ====== PARIS GDML (inchangé) ======
    G4GDMLParser parser;
    parser.Read("../PARISMalia.gdml");
    auto lvHousing  = parser.GetVolume("SCIONIXPWLVFullHousing");
    auto lvCe       = parser.GetVolume("SCIONIXPWLVCe");
    auto lvNaI      = parser.GetVolume("SCParisPWLV.1");
    auto lvQuartzPM = parser.GetVolume("SCParisPWQuartzforPMLV");
    auto lvSealPM   = parser.GetVolume("SCParisPWSealingforPMLV");
    if(!(lvHousing && lvCe && lvNaI && lvQuartzPM && lvSealPM)) {
        G4Exception("MyDetectorConstruction::Construct","PARIS-MissingLV",FatalException,
                    "Logical(s) PARIS introuvable(s) dans le GDML.");
    }
    auto asmPARIS = new G4AssemblyVolume();
    G4RotationMatrix rotId;
    G4ThreeVector trHousing ( 193.*mm, 0., 230.5*mm );
    G4ThreeVector trCe      ( 193.*mm, 0., 258.5*mm );
    G4ThreeVector trNaI     ( 193.*mm, 0., 360.01*mm );
    G4ThreeVector trQuartz  ( 193.*mm, 0., 439.02*mm );
    G4ThreeVector trSeal    ( 193.*mm, 0., 436.02*mm );
    G4Transform3D T_housing ( rotId, trHousing );
    G4Transform3D T_ce      ( rotId, trCe      );
    G4Transform3D T_nai     ( rotId, trNaI     );
    G4Transform3D T_quartz  ( rotId, trQuartz  );
    G4Transform3D T_seal    ( rotId, trSeal    );
    asmPARIS->AddPlacedVolume(lvHousing,  T_housing);
    asmPARIS->AddPlacedVolume(lvCe,       T_ce);
    asmPARIS->AddPlacedVolume(lvNaI,      T_nai);
    asmPARIS->AddPlacedVolume(lvQuartzPM, T_quartz);
    asmPARIS->AddPlacedVolume(lvSealPM,   T_seal);

    //const G4ThreeVector arcCenter(-60.5*4*mm, 0., 0.);
    const G4ThreeVector arcCenter(-40*4-33*mm, 0., 0.);
    const std::array<G4double,9> thetas = {50.*deg, 70.*deg, 90.*deg, 110.*deg, 130.*deg, 235.*deg, 262.*deg, 278.*deg, 305.*deg};
    int copyNo = 0;
    const G4RotationMatrix RA_top = *rotY;
    const G4ThreeVector Cworld(0.,0.,0.);
    const G4double D_frontCe = 233.0*mm; //+6.656mm peut-être à ajuster selon le GDML
    //const G4double D_frontCe = 207.5*mm; // test d'après le GDML fourni
    const G4double Rtarget   = 300.0*mm;
    for (auto theta : thetas) {
        const G4double phi_loc = theta;
        G4ThreeVector r_local(std::cos(phi_loc), std::sin(phi_loc), 0.);
        G4ThreeVector w = RA_top * r_local;  w = w.unit();
        const G4ThreeVector pos = Cworld + ((Rtarget - D_frontCe) * w) + arcCenter;
        //const G4ThreeVector pos = Cworld + (Rtarget + D_frontCe) * w;        // test(sans arcCenter)
        //const G4ThreeVector pos = Cworld - arcCenter + (Rtarget + D_frontCe) * w;
        const G4ThreeVector k(0.,0.,1.);
        G4double cang = std::max(-1.0, std::min(1.0, (double)k.dot(w)));
        G4double ang = std::acos(cang);
        G4ThreeVector axis = k.cross(w);
        if (axis.mag2() < 1e-24) axis = G4ThreeVector(1,0,0);
        axis = axis.unit();
        G4RotationMatrix R;  R.rotate(ang, axis);
        G4Transform3D T(R, pos);
        asmPARIS->MakeImprint(logicWorld, T, copyNo, checkOverlaps);

        const int degLab = (int) std::round(theta/deg);
        SetParisLabel(copyNo,"PARIS" + std::to_string(degLab)); // ex: PARIS50, PARIS90, ...
        copyNo++;

        G4ThreeVector zplus_world = R * k;
        G4ThreeVector faceCe = pos - D_frontCe * zplus_world;
        G4double dist = faceCe.mag();
        G4cout << "θ=" << theta/deg << "°  dist(faceCe->world) = " << dist/mm << " mm" << G4endl;
    }

    // ====== Retour ======
    return physWorld;
}

void MyDetectorConstruction::ConstructSDandField()
{
    auto sdMan = G4SDManager::GetSDMpointer();

    // Ce & NaI (énergie déposée)
    auto lvCe  = G4LogicalVolumeStore::GetInstance()->GetVolume("SCIONIXPWLVCe");
    auto lvNaI = G4LogicalVolumeStore::GetInstance()->GetVolume("SCParisPWLV.1");

    auto sdCe = new G4MultiFunctionalDetector("CeSD");
    sdMan->AddNewDetector(sdCe);
    // NOTE: depth=0 pour indexer par le copy number du parent (imprint PARIS)
    auto edepCe = new G4PSEnergyDeposit("edep", /*depth=*/0);
    sdCe->RegisterPrimitive(edepCe);
    if (lvCe) lvCe->SetSensitiveDetector(sdCe);

    auto sdNaI = new G4MultiFunctionalDetector("NaISD");
    sdMan->AddNewDetector(sdNaI);
    // NOTE: depth=0 pour indexer par le copy number du parent (imprint PARIS)
    auto edepNaI = new G4PSEnergyDeposit("edep", /*depth=*/0);
    sdNaI->RegisterPrimitive(edepNaI);
    if (lvNaI) lvNaI->SetSensitiveDetector(sdNaI);

    // Cells : filtre neutron E < 100 keV
    // auto filterThermalN = new G4SDParticleWithEnergyFilter("fThermalN");
    // filterThermalN->add("neutron");
    // filterThermalN->SetKineticEnergy(0.*eV, 10000.*keV);

    auto sdCell = new G4MultiFunctionalDetector("CellSD");
    sdMan->AddNewDetector(sdCell);

    auto psEnter = new G4PSTrackCounter("nThermalEnter", fCurrent_In);
    //psEnter->SetFilter(filterThermalN);
    sdCell->RegisterPrimitive(psEnter);

    if (fScoringVolumeOne)   fScoringVolumeOne->SetSensitiveDetector(sdCell);
    if (fScoringVolumeTwo)   fScoringVolumeTwo->SetSensitiveDetector(sdCell);
    if (fScoringVolumeThree) fScoringVolumeThree->SetSensitiveDetector(sdCell);
    if (fScoringVolumeFour)  fScoringVolumeFour->SetSensitiveDetector(sdCell);
}