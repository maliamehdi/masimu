#include "DetectorConstruction.hh"

void PlaceRingCells(
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
    
    Pressure = 7*6.24151e+08; // MeV/mm3
}

MyDetectorConstruction::~MyDetectorConstruction()
{}

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{	
    auto nist = G4NistManager::Instance();

    // Matériaux
    auto air = nist->FindOrBuildMaterial("G4_AIR");
    auto poly = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
    auto boron = nist->FindOrBuildMaterial("G4_B");
    auto methane = nist->FindOrBuildMaterial("G4_METHANE");
    // Définition ABS (Acrylonitrile Butadiène Styrène) Plastique impression 3D
    G4double density = 1.04*g/cm3;
    G4int ncomponents;
    G4Material* ABS = new G4Material("ABS", density, ncomponents=3);

    G4Element* elC = new G4Element("Carbon", "C", 6., 12.011*g/mole);
    G4Element* elH = new G4Element("Hydrogen", "H", 1., 1.008*g/mole);
    G4Element* elN = new G4Element("Nitrogen","N", 7., 14.007*g/mole);

    ABS->AddElement(elC, 85*perCent);
    ABS->AddElement(elH, 9*perCent);
    ABS->AddElement(elN, 6*perCent);
    //Steel cylinder
    steel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
   //Helium 3
	G4double atomicMass = 3.016*g/mole;
	
	G4Isotope* he3Iso = new G4Isotope("he3Iso", 2, 3, atomicMass);
	
	G4Element *he3 = new G4Element("he3", "he3", 1);
	
	he3->AddIsotope(he3Iso, 100.*perCent);

    //Cell gas mixture
    G4double    temperatureCell = 293.*kelvin,
                molar_constant = CLHEP::Avogadro*CLHEP::k_Boltzmann;
	
	G4double    pressureCellOne = 2*bar,//4.8*atmosphere, //mettre à 2bar ça marche du feu de dieu pour Emile
                densityCellOne = (atomicMass*pressureCellOne)/(temperatureCell*molar_constant);

    G4double    pressureCellTwo = 2*bar,//5.6*atmosphere,
                densityCellTwo = (atomicMass*pressureCellTwo)/(temperatureCell*molar_constant);

    G4double    pressureCellThree = 2*bar,// 5.8*atmosphere,
                densityCellThree = (atomicMass*pressureCellThree)/(temperatureCell*molar_constant);

    G4double    pressureCellFour = 2*bar,//5.8*atmosphere,
                densityCellFour = (atomicMass*pressureCellFour)/(temperatureCell*molar_constant);
                
    G4Material* co2 = nist->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
		        
	G4Material *gasMixOne = new G4Material("he3GasOne", densityCellOne,	2, kStateGas, temperatureCell, pressureCellOne);
    
    gasMixOne->AddElement(he3, 99.*perCent);
    gasMixOne->AddMaterial(co2, 1.*perCent);
    
    G4Material* gasMixTwo = new G4Material("he3GasTwo", densityCellTwo,	2, kStateGas, temperatureCell, pressureCellTwo);
    
    gasMixTwo->AddElement(he3, 99.*perCent);
    gasMixTwo->AddMaterial(co2, 1.*perCent);

    G4Material* gasMixThree = new G4Material("he3GasThree", densityCellThree,	2, kStateGas, temperatureCell, pressureCellThree);
    
    gasMixThree->AddElement(he3, 99.*perCent);
    gasMixThree->AddMaterial(co2, 1.*perCent);

    G4Material* gasMixFour = new G4Material("he3GasFour", densityCellFour,	2, kStateGas, temperatureCell, pressureCellFour);
    
    gasMixFour->AddElement(he3, 99.*perCent);
    gasMixFour->AddMaterial(co2, 1.*perCent);
	    //Aluminium
        G4Material* aluminium = nist->FindOrBuildMaterial("G4_Al");


        modMat = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
        //High density polyethylene (moderator)
        mod = new G4Material("mod", 0.93*g/cm3, 1);
        mod->AddMaterial(modMat, 100.*perCent);
         //Borated polyethylène
        borPol = new G4Material("borPol", 1.005*g/cm3, 2);
        borPol->AddMaterial(poly, 95.*perCent);
        borPol->AddMaterial(boron, 5.*perCent);
    
        //Plastic beta detector
        plasticMat = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    
        //Steel cylinder
        steel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    
        //Plexiglass
        plexi = nist->FindOrBuildMaterial("G4_PLEXIGLASS");


        // 1. Boite d'air
        solidWorld = new G4Box("solidWorld", 1.*m, 1.*m, 1.*m);
        logicWorld = new G4LogicalVolume(solidWorld, air, "logicWorld");
        physWorld = new G4PVPlacement(0, {}, logicWorld, "physWorld", nullptr, false, 0);
        
        //========== TETRA MODERATOR ==========
        // Parameters
        const G4double R = 36.6*cm;      // Hexagon "radius" core polyethylène
        const G4double Rshield = 51.6*cm; // Hexagon "radius" shielding shell borated polyethylene
        const G4double halfHeight = 25.*cm;      // Half-height of extrusion
        const G4double offset = 148.5*mm;        // Half-gap between halves
        const G4double holeRadius = 6.75*cm;     // Hole radius


        //---------------------------------------------------
        // =======A laMatthieu avec Polyhedra========

        G4bool checkOverlaps = false;
        // Demi-hexagone asymétrique
        // I define The first half of hexagone
        // From front view:
        // Left volumes
        G4int numRZ = 4;
        G4double rpos[]={0*mm,516.*mm,516*mm,0*mm};//,366.*mm,19.*mm,19.*mm,366.*mm,367.*mm,364.*mm,363.*mm,367.*mm,368.*mm,367.*mm};
        G4double zpos[]={-350*mm,-350.*mm,350.*mm,350*mm};
        G4double rpos2[]={0*mm,366.*mm,366*mm,0*mm};
        G4double zpos2[]={-250.1*mm,-250.1*mm,250.1*mm,250.1*mm};
        // Volume for borated polypropylene
        G4Polyhedra *solid_half_shell =   new G4Polyhedra( "solid_shell",
                                                    270*deg,           //- initial phi starting angle
                                                    180.*deg,    //- total phi angle
                                                    3,           //- number sides : 9 on sides + 2 front and back
                                                    numRZ,          //- number corners in r,z space
                                                    rpos,        //- r coordinate of these corners
                                                    zpos );      //- z coordinate of these corners
        G4Polyhedra *solid_half_shell2 = new G4Polyhedra("solid_shell2",
                                                  90*deg,         // φstart symétrique
                                                  180.*deg,       // φtotal
                                                  3,
                                                  numRZ,
                                                  rpos,
                                                  zpos);
        // Getting the same hexagone but smaller
        G4Polyhedra* solid_half_polycase  = new G4Polyhedra("solid_half_polycase", 270*deg, 180.*deg, 3, numRZ, rpos2, zpos2);
        G4Polyhedra *solid_half_polycase2 = new G4Polyhedra("solid_polycase2",
                                                     90*deg,        // φstart inversé
                                                     180.*deg,
                                                     3,
                                                     numRZ,
                                                     rpos2,
                                                     zpos2);
        // Add-ons (boîtes pour créer la forme asymétrique)
        G4Box* box_addon_shell     = new G4Box("shell_addon", 50./3.*mm, 258*mm, 350*mm);
        G4Box* box_addon_polycase  = new G4Box("polycase_addon", 50./3.*mm, (183-16.875)*mm, 250*mm);
        
        // Union + Subtraction : shell 1
        G4UnionSolid* solid_shell_init1 = new G4UnionSolid("solid_shell_init1", solid_half_shell, box_addon_shell, nullptr, G4ThreeVector(-50./3./2.*mm, -258*mm, 0));
        G4SubtractionSolid* solid_shell_init = new G4SubtractionSolid("solid_shell_init", solid_shell_init1, box_addon_shell, nullptr, G4ThreeVector(-50./3./2.*mm, 258*mm, 0));

        // Union + Subtraction : polycase 1
        G4UnionSolid* solid_polycase_init1 = new G4UnionSolid("solid_polycase_init1", solid_half_polycase, box_addon_polycase, nullptr, G4ThreeVector(-50./3./2.*mm, (-(183+16.875))*mm, 0));
        G4SubtractionSolid* solid_polycase_init = new G4SubtractionSolid("solid_polycase_init", solid_polycase_init1, box_addon_polycase, nullptr, G4ThreeVector(-50./3./2.*mm, (183+16.875)*mm, 0));

        // Union + Subtraction : shell 2 (symétrique)
        G4UnionSolid* solid_shell2_init1 = new G4UnionSolid("solid_shell2_init1", solid_half_shell2, box_addon_shell, nullptr, G4ThreeVector(+50./3./2.*mm, 258*mm, 0));
        G4SubtractionSolid* solid_shell2_init = new G4SubtractionSolid("solid_shell2_init", solid_shell2_init1, box_addon_shell, nullptr, G4ThreeVector(+50./3./2.*mm, -258*mm, 0));

        // Union + Subtraction : polycase 2 (symétrique)
        G4UnionSolid* solid_polycase2_init1 = new G4UnionSolid("solid_polycase2_init1", solid_half_polycase2, box_addon_polycase, nullptr, G4ThreeVector(+50./3./2.*mm, (183+16.875)*mm, 0));
        //G4SubtractionSolid* solid_polycase2_init2 = new G4SubtractionSolid("solid_polycase2_init2", solid_polycase2_init1, box_cut, nullptr, G4ThreeVector(+50./3./2.*mm, 0, 0));
        G4SubtractionSolid* solid_polycase2_init = new G4SubtractionSolid("solid_polycase2_init", solid_polycase2_init1, box_addon_polycase, nullptr, G4ThreeVector(+50./3./2.*mm, (-(183+16.875))*mm, 0));

        // Soustraction du tuyau central
        G4VSolid* solid_hole = new G4Tubs("solid_pipe", 0*mm, 67.5*mm, 520.*mm, 0.*deg, 360.0*deg);
        // Creating the light guide cut
        G4VSolid* solidLightGuide = new G4Box("solidLightGuide", 2.25*cm, 25.*cm, 5.*cm);
        G4ThreeVector yTrans(0., -20.*cm, 0.);
        // Subtracting the hole tube the light guide from the shell and polycase
        G4UnionSolid* solid_pipe = new G4UnionSolid("solid_pipe", solid_hole, solidLightGuide, nullptr, yTrans);
        G4SubtractionSolid* solid_shell = new G4SubtractionSolid("solid_shell", solid_shell_init, solid_pipe);
        G4SubtractionSolid* solid_polycase = new G4SubtractionSolid("solid_polycase", solid_polycase_init, solid_pipe);
        G4SubtractionSolid* solid_shell2 = new G4SubtractionSolid("solid_shell2", solid_shell2_init, solid_pipe);
        G4SubtractionSolid* solid_polycase2 = new G4SubtractionSolid("solid_polycase2", solid_polycase2_init, solid_pipe);

        // Volumes logiques
        G4LogicalVolume* logic_shell = new G4LogicalVolume(solid_shell, borPol, "logic_shell");
        G4LogicalVolume* logic_polycase = new G4LogicalVolume(solid_polycase, mod, "logic_polycase");
        G4LogicalVolume* logic_shell2 = new G4LogicalVolume(solid_shell2, borPol, "logic_shell2");
        G4LogicalVolume* logic_polycase2 = new G4LogicalVolume(solid_polycase2, mod, "logic_polycase2");

        
        
        G4LogicalVolume *logicLightGuide = new G4LogicalVolume(solidLightGuide,
                                                                nist->FindOrBuildMaterial("air"),
                                                                "logicLightGuide");
        //To be placed later in the final logical volumes
        G4ThreeVector origin = G4ThreeVector(0,0,0);


        //======Partie pour les cellules de gaz========
        //G4ThreeVector origin = G4ThreeVector(0,0,0);
        // Now I create two tubes:
        // One of steel as the gas cell
        // and another one made of gas
        //Cell cylinder
	    G4VSolid* solidCyl = new G4Tubs("CylOut", 15.5*mm, 16.*mm, 250.*mm, 0., 360.);
	    G4LogicalVolume* logicCyl = new G4LogicalVolume(solidCyl, steel, "logicCyl");
        //Helium cells
        G4VSolid* solidCell = new G4Tubs("solidCell", 0., 15.5*mm, 250.*mm, 0., 360.);
        
        G4LogicalVolume* logicCellOne = new G4LogicalVolume(solidCell, gasMixOne, "logicCellOne");
        G4LogicalVolume* logicCellTwo = new G4LogicalVolume(solidCell, gasMixTwo, "logicCellTwo");
        G4LogicalVolume* logicCellThree = new G4LogicalVolume(solidCell, gasMixThree, "logicCellThree");
        G4LogicalVolume* logicCellFour = new G4LogicalVolume(solidCell, gasMixFour, "logicCellFour");
        
        fScoringVolumeOne = logicCellOne;
        fScoringVolumeTwo = logicCellTwo;
        fScoringVolumeThree = logicCellThree;
        fScoringVolumeFour = logicCellFour;
        
        //Cells physical volume
        G4int nsect1 = 6;
        G4double dphi1 = 360.*deg;
        G4double ang1 = dphi1/nsect1;
        std::vector<G4ThreeVector> hex(nsect1);
    
        G4double c = 5.*cm; // comme avant
        G4int copyIndex = 0;
        PlaceRingCells(2.*c,      -30.*deg, 6, logicCellOne,   logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
        PlaceRingCells(sqrt(3.)*c,  0.*deg, 6, logicCellOne,   logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
        PlaceRingCells(3.*c,      -30.*deg, 6, logicCellTwo,   logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
        PlaceRingCells(sqrt(7.)*c, 10.893*deg, 6, logicCellTwo, logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
        PlaceRingCells(sqrt(7.)*c,-10.893*deg, 6, logicCellTwo, logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
        PlaceRingCells(4.*c,      -30.*deg, 6, logicCellThree, logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
        PlaceRingCells(2.*sqrt(3.)*c,  0.*deg, 6, logicCellThree, logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
        PlaceRingCells(sqrt(13.)*c, 16.102*deg, 6, logicCellThree, logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
        PlaceRingCells(sqrt(13.)*c,-16.102*deg, 6, logicCellThree, logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
        PlaceRingCells(5.*c,      -30.*deg, 6, logicCellFour,  logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
        PlaceRingCells(sqrt(19.)*c, 6.587*deg, 6, logicCellFour,  logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
        PlaceRingCells(sqrt(21.)*c, 19.107*deg, 6, logicCellFour,  logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
        PlaceRingCells(sqrt(21.)*c,-19.107*deg, 6, logicCellFour,  logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
        PlaceRingCells(sqrt(19.)*c,-6.587*deg, 6, logicCellFour,  logicCyl, logic_polycase, logic_polycase2, origin, origin, copyIndex);
        
            
        //--------------------Fin placement des tubes--------------------
        // Placement des volumes physiques chargés de cellules de gaz
        G4ThreeVector shell_pos = G4ThreeVector(+offset,0,0); //to open TETRA
        G4ThreeVector shell2_pos = G4ThreeVector(-offset,0,0);
        G4VPhysicalVolume *phys_polycase = new G4PVPlacement(0,
                                 origin,
                                 logic_polycase,
                                 "Polycase",
                                 logic_shell,
                                 false,0,checkOverlaps);
        G4VPhysicalVolume *phys_polycase2 = new G4PVPlacement(0,
                                 origin,
                                 logic_polycase2,
                                 "Polycase",
                                 logic_shell2,
                                 false,0,checkOverlaps);

        G4VPhysicalVolume *phys_shell = new G4PVPlacement(0,
                                        shell_pos,
                                        logic_shell,
                                        "Shell",
                                        logicWorld,
                                        false,0,checkOverlaps);
        G4VPhysicalVolume *phys_shell2 = new G4PVPlacement(0,
                                        shell2_pos,
                                        logic_shell2,
                                        "Shell2",
                                        logicWorld,
                                        false,0,checkOverlaps);
       
        //---------------------------------------------------
        //Twin Frisch-grid ionization chamber
        // 1. Create the solid for the ionization chamber
        //G4Tubs* solidIC = new G4Tubs("solidIC", 0., 52.*mm, 198.5*mm, 0., 360);
        // G4VSolid* solidVac = new G4Tubs("solidVac", 0., 52.*mm, 198.5*mm, 0., 360.);
        
        // G4LogicalVolume* logicVac = new G4LogicalVolume(solidVac, air, "logicVac");
        
        // new G4PVPlacement(0, G4ThreeVector(0., 0., -151.5*mm), logicVac, "physVac", logicUni, false, 0, true);

    
        //---------------------------------------------------
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

        // base = ancien offset du mother : (0, +590, 0)
        const G4ThreeVector base(0., 590.*mm, 0.);

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
        new G4PVPlacement(nullptr, base + G4ThreeVector(0,  -838*mm, -130.6*mm),
                        lv_barre_5, "pv_barre_5_L", logicWorld, false, 0, checkOverlaps);
        new G4PVPlacement(nullptr, base + G4ThreeVector(0,  -838*mm, +138.9*mm),
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
        G4PVPlacement* physArc1 = new G4PVPlacement(rotY, G4ThreeVector(-60.5*mm, 0.*mm, 0.*mm), logicArc, "physArc1", logicWorld, false, 0, false);
        G4PVPlacement* physArc2 = new G4PVPlacement(rotY2, G4ThreeVector(-60.5*mm, 0.*mm, 0.*mm), logicArc, "physArc2", logicWorld, false, 1, false);

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
        const G4double w_arc_supp    = 30.*mm;   // EXEMPLE: 30 mm (met ta valeur)

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

        // --- Variante : si tu veux l’axe du cylindre = Y (au lieu de Z) ---
        // G4RotationMatrix* rotX90 = new G4RotationMatrix(); rotX90->rotateX(90.*deg);
        // new G4PVPlacement(rotX90, G4ThreeVector(0,+50*mm,0),
        //                   logicSupportWithArc, "physSupportABS_A",
        //                   logicWorld, false, 0, true);
        // new G4PVPlacement(rotX90, G4ThreeVector(0,-50*mm,0),
        //                   logicSupportWithArc, "physSupportABS_B",
        //                   logicWorld, false, 1, true);
        //---------------------------------------------------
        //Chambre d'ionisation
        G4double rICint = 89.5*mm; // Rayon interne de la chambre d'ionisation
        G4double rICext = 90.*mm; // Rayon externe de la chambre d'ionisation
        G4double hIC = 77.675*mm; // Demi hauteur de la chambre d'ionisation

        G4Tubs* solidIC = new G4Tubs("solidIC", rICint, rICext, hIC, 0., 360.*deg);
        //rajouter un cylindre plein solide au bout de la chamnbre d'ionisation
        G4Tubs* solidICend = new G4Tubs("solidICend", 0., rICext+10*mm, 13.37*mm, 0., 360.*deg);
        //Unir le cylindre plein au bout de la chambre d'ionisation
        G4UnionSolid* solidICfinal = new G4UnionSolid("solidICfinal", solidIC, solidICend, nullptr, G4ThreeVector(0., 0., hIC + 13.37*mm));
        //Créer le volume logique de la chambre d
        G4LogicalVolume* logicIC = new G4LogicalVolume(solidICfinal, aluminium, "logicIC");
        // Cylindre de Methane à l'intérieur de la chambre d'ionisation
        G4Tubs* solidMethane = new G4Tubs("solidMethane", 0., rICint, hIC, 0., 360.*deg);
        G4LogicalVolume* logicMethane = new G4LogicalVolume(solidMethane, methane, "logicMethane");
        //Placer le cylindre de methane dans la chambre d'ionisation
        G4VPhysicalVolume* physMethane = new G4PVPlacement(0,
                                        G4ThreeVector(0., 0., 0.),
                                        logicMethane,
                                        "physMethane",
                                        logicIC,
                                        false, 0, checkOverlaps);
        //Placer la chambre d'ionisation dans le monde
        G4VPhysicalVolume* physIC = new G4PVPlacement(0,
                                        origin,
                                        logicIC,
                                        "physIC",
                                        logicWorld,
                                        false, 0, checkOverlaps);

        

        //========Source Point========
        //Source Point
 	    solidSP = new G4Sphere("SP", 0., 0.1*mm, 0., 360., 0., 180.);
 	    logicSP = new G4LogicalVolume(solidSP, air, "logicSP");
 	    physSP = new G4PVPlacement(0, G4ThreeVector(0., -2.*cm, 10.*cm), logicSP, "physSP", logicWorld, false, 0, true);
	
 	    fSPVolume = logicSP;
        
        //==========PARIS CeBr3/NaI DETECTORS==========
        // Constructing PARIS phoswiches //
        // === Import du GDML ===
            // ========== PARIS CeBr3/NaI (via GDML) – Assembly + placement sur l’arc haut ==========
        G4GDMLParser parser;
        parser.Read("../PARISMalia.gdml");

        // Logical volumes composant un phoswich PARIS
        auto lvHousing  = parser.GetVolume("SCIONIXPWLVFullHousing");
        auto lvCe       = parser.GetVolume("SCIONIXPWLVCe");
        auto lvNaI      = parser.GetVolume("SCParisPWLV.1");
        auto lvQuartzPM = parser.GetVolume("SCParisPWQuartzforPMLV");
        auto lvSealPM   = parser.GetVolume("SCParisPWSealingforPMLV");

        if(!(lvHousing && lvCe && lvNaI && lvQuartzPM && lvSealPM)) {
        G4Exception("MyDetectorConstruction::Construct","PARIS-MissingLV",FatalException,
                    "Logical(s) PARIS introuvable(s) dans le GDML.");
        }

        // Assembly = “recette” d’un module
        auto asmPARIS = new G4AssemblyVolume();

        // ---- Positions RELATIVES (bloc .000 du GDML) ----
        // Option A: via Transform3D **lvalues**
        G4RotationMatrix rotId; // identité
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

        // --- Placement sur l’arc du haut ---
        // Repère : l’arc (logicArc) est dans le plan (Y,Z). Angle θ mesuré de +Y vers +Z.
        // Vecteur radial sortant : w = (0, cosθ, sinθ).
        // On aligne +Z(local) sur w ⇒ la face Ce (−Z local) regarde le centre.
        const G4ThreeVector arcCenter(-60.5*4*mm, 0., 0.);  // centre de ton arc haut
        
        const std::array<G4double,5> angles = {50.*deg, 70.*deg, 90.*deg, 110.*deg, 130.*deg};
        // --- même rotation que le rail du haut ---
        const G4RotationMatrix RA_top = *rotY;          // rot du rail (Y=+90°)
        const G4ThreeVector    C(0.,0.,0.);             // centre = world
        const G4double         D_frontCe = 233.0*mm;    // origine module -> face Ce (côté centre)
        const G4double         Rtarget   = 300.0*mm;    // distance voulue face Ce -> centre world

        const std::array<G4double,9> thetas = {50.*deg, 70.*deg, 90.*deg, 110.*deg, 130.*deg, 235.*deg, 262.*deg, 278.*deg, 305.*deg};
        int copyNo = 0;

        for (auto theta : thetas) {
            // angle local du Tubs après rotY=+90° : φ_local = 90° + θ
            const G4double phi_loc =  theta;

            // radial LOCAL (dans le plan XY natif du Tubs)
            G4ThreeVector r_local(std::cos(phi_loc), std::sin(phi_loc), 0.);
            // radial MONDE (après rotation du rail) et normalisé
            G4ThreeVector w = RA_top * r_local;  w = w.unit();

            // position de l’ORIGINE du module pour que la face Ce soit à 300 mm du centre world
            const G4ThreeVector pos = C + ((Rtarget - D_frontCe) * w) + arcCenter; 

            // orienter +Z(local) du module sur w  (=> la face Ce, côté -Z, regarde le centre)
            const G4ThreeVector k(0.,0.,1.);
            G4double c = std::max(-1.0, std::min(1.0, (double)k.dot(w)));
            G4double ang = std::acos(c);
            G4ThreeVector axis = k.cross(w);
            if (axis.mag2() < 1e-24) axis = G4ThreeVector(1,0,0);
            axis = axis.unit();

            G4RotationMatrix R;  R.rotate(ang, axis);

            // Imprint (API 11.2.x : passer des lvalues)
            G4Transform3D T(R, pos);
            asmPARIS->MakeImprint(logicWorld, T, copyNo++, checkOverlaps);

            // --- Vérif : distance face Ce -> centre world ~ 300 mm ---
            G4ThreeVector zplus_world = R * k;                // +Z local dans le monde
            G4ThreeVector faceCe = pos - D_frontCe * zplus_world;  // côté -Z
            G4double dist = faceCe.mag();                     // centre world = (0,0,0)
            G4cout << "θ=" << theta/deg << "°  dist(faceCe->world) = "
                    << dist/mm << " mm" << G4endl;
        } // for angles
        // === Fin de l'import du GDML ===
       
       
   
    
            //Mono Cell for tests
    /*monoCell = new G4PVPlacement(0, G4ThreeVector(0., (2.*c)*cm, 0.), logicCell, "physCell", logicMod, false, 0, true);
    monoCyl = new G4PVPlacement(0, G4ThreeVector(0., (2.*c)*cm, 0.), logicCyl, "physCyl", logicMod, false, 0, true);*/
// Exemple simplifié de détecteur composite hexagonal coupé en deux moitiés (ModA et ModB)
// avec espace de 29.7 cm dans une boite d'air
return physWorld;
}