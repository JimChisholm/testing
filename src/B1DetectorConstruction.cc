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
// $Id: B1DetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //
  G4double env_sizeXY = 262*cm, env_sizeZ = 262*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Envelope
  //  
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 
  //     
  // Shape 1
  //  
  //G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_CONCRETE");
 // G4ThreeVector pos1 = G4ThreeVector(1, 3*cm, -11*cm);       
  // Conical section shape       
  //G4double shape1_rmina =  0.4*cm, shape1_rmaxa = 1.2*cm;
 // G4double shape1_rminb =  0.3*cm, shape1_rmaxb = 1.4*cm;
  //G4double shape1_hz = 2.*cm;
  //G4double shape1_phimin = 150.*deg, shape1_phimax = 260.*deg;
 // G4Cons* solidShape1 =    
    //new G4Cons("Shape1", 
    //shape1_rmina, shape1_rmaxa, shape1_rminb, shape1_rmaxb, shape1_hz,
    //shape1_phimin, shape1_phimax);                      
 // G4LogicalVolume* logicShape1 =                         
    //new G4LogicalVolume(solidShape1,         //its solid
          //              shape1_mat,          //its material
           //             "Shape1");           //its name
               
  //new G4PVPlacement(0,                       //no rotation
                    //pos1,                    //at position
                    //logicShape1,             //its logical volume
                    //"Shape1",                //its name
                    //logicEnv,                //its mother  volume
                   // false,                   //no boolean operation
                    //0,                       //copy number
                   // checkOverlaps);          //overlaps checking

G4Element* Oxygen = new G4Element("Oxygen","O", 8, 16*g/mole);
G4Element* Sulfur = new G4Element("Sulfur","S", 16, 32.06*g/mole);
G4Element* Silicon = new G4Element("Silicon","Si", 14, 28.085*g/mole);
G4Element* Aluminum = new G4Element("Aluminum","Al", 13, 26.982*g/mole);
G4Element* Iron = new G4Element("Iron","Fe", 26, 55.845*g/mole);
G4Element* Magnesium = new G4Element("Magnesium","Mg", 12, 24.305*g/mole);
G4Element* Manganese = new G4Element("Manganese", "Mn", 25, 54.938*g/mole);
G4Element* Calcium = new G4Element("Calcium","Ca", 20, 40.078*g/mole);
G4Element* Sodium = new G4Element("Sodium","Na", 11, 22.990*g/mole);
G4Element* Potassium = new G4Element("Potassium","K", 19, 26.982*g/mole);
G4Element* Phosphorus = new G4Element("Phosphorus","P", 15, 30.974*g/mole);
G4Element* Chlorine = new G4Element("Chlorine","Cl", 17, 35.45*g/mole);
//G4Element* Nitrogen = new G4Element("Nitrogen","N", 7, 14.007*g/mole);
//G4Element* Zinc = new G4Element("Zinc","Zn", 30, 65.38*g/mole);
G4Element* Chromium = new G4Element("Chromium","Cr", 17, 35.45 *g/mole);
G4Element* Titanium = new G4Element("Titanium","Ti", 22, 44.956*g/mole);
G4Element* Hydrogen = new G4Element("Hydrogen", "H", 1, 1.008*g/mole);
G4Element* Carbon = new G4Element("Carbon", "C", 6, 12.011*g/mole);

//G4Material* CL = new G4Material("Chlorine", 3.2*g/cm3, 2); //do not know densities for sure as these are stp.
//CL->AddElement(Chlorine, 1);

G4Material* SIO2 = new G4Material("Silicondioxide", 2.195*g/cm3, 2); 
SIO2->AddElement(Silicon, 1);
SIO2->AddElement(Oxygen, 2);

G4Material* Cr2O3 = new G4Material("Chromiumtrioxide", 5.22*g/cm3, 2);
Cr2O3->AddElement(Chromium, 2);
Cr2O3->AddElement(Oxygen, 3);

G4Material* P2O5 = new G4Material("PhosphorusPentoxide", 2.39*g/cm3, 2);
P2O5->AddElement(Phosphorus, 2);
P2O5->AddElement(Oxygen, 5);

G4Material* K2O = new G4Material("PotassiumOxide", 2.32*g/cm3, 2);
K2O->AddElement(Potassium, 2);
K2O->AddElement(Oxygen, 1);

G4Material* Na2O = new G4Material("SodiumDioxide", 5.745*g/cm3, 2);
Na2O->AddElement(Sodium, 2);
Na2O->AddElement(Oxygen, 1);

G4Material* CaO = new G4Material("QuickLime", 3.34*g/cm3, 2);
CaO->AddElement(Calcium, 1);
CaO->AddElement(Oxygen, 1);

G4Material* MnO = new G4Material("ManganeseOxide", 5.37*g/cm3, 2);
MnO->AddElement(Manganese, 1);
MnO->AddElement(Oxygen, 1);

G4Material* SO3 = new G4Material("SulfurTrioxide", 1.92*g/cm3, 2); 
SO3->AddElement(Sulfur, 1);
SO3->AddElement(Oxygen, 3);

G4Material* TiO2 = new G4Material("TitaniumDioxide", 4.23*g/cm3, 2);
TiO2->AddElement(Titanium, 1);
TiO2->AddElement(Oxygen, 2);

G4Material* FeO = new G4Material("IronOxide", 5.745*g/cm3, 2);
FeO->AddElement(Iron, 1);
FeO->AddElement(Oxygen, 1);

G4Material* Fe2O3 = new G4Material("Fe3Oxide", 5.25*g/cm3, 2);
Fe2O3->AddElement(Iron, 2);
Fe2O3->AddElement(Oxygen, 3);

G4Material* MgO = new G4Material("MagnesiumDioxide", 5.43*g/cm3, 2);
MgO->AddElement(Magnesium, 1);
MgO->AddElement(Oxygen, 1);

G4Material* Al2O3 = new G4Material("AluminumOxide", 3.95*g/cm3, 2);
Al2O3->AddElement(Aluminum, 2);
Al2O3->AddElement(Oxygen, 3);


//shape3
//G4Material* shape3_mat = nist->FindOrBuildMaterial("G4_CONCRETE");
//G4Material* Al2CaK2MgNa2O7 = new G4Material("MarsRock", 1.62*g/cm3, 6); //needs density for Al2CaK2MgNa2O7
//Al2CaK2MgNa2O7->AddElement(Aluminum, 2);
//Al2CaK2MgNa2O7->AddElement(Calcium, 1);
//Al2CaK2MgNa2O7->AddElement(Potassium, 2);
//Al2CaK2MgNa2O7->AddElement(Magnesium, 1);
//Al2CaK2MgNa2O7->AddElement(Sodium, 2);
//Al2CaK2MgNa2O7->AddElement(Oxygen, 7);

//Polyimide H: 1.58 atom density 10^22atoms/gm; O: 0.788 atomic density. 1.42 density
//Polyimide C 3.47 atomic density; N 0.315 atomic density Kim 1994 p. 14 

G4Material*MMS2 = new G4Material("MarsRock", 1.62*g/cm3, 12);
MMS2->AddMaterial(SIO2, 43.8*perCent);
//shape2_mat->AddMaterial(FeO, 10.*perCent);
MMS2->AddMaterial(Fe2O3, 18.37*perCent); //Used internet source bc sara's is not complete
MMS2->AddMaterial(MnO, 0.13*perCent); 
MMS2->AddMaterial(MgO, 6.66*perCent);
MMS2->AddMaterial(CaO, 7.98*perCent);
MMS2->AddMaterial(Na2O, 2.51*perCent);
MMS2->AddMaterial(K2O, 0.37*perCent);
MMS2->AddMaterial(P2O5, 0.13*perCent);
MMS2->AddMaterial(Cr2O3, 0.04*perCent);
MMS2->AddMaterial(SO3, 6.11*perCent);
MMS2->AddMaterial(TiO2, 0.83*perCent);
MMS2->AddMaterial(Al2O3, 13.07*perCent);
//shape2_mat->AddMaterial(CL, 10.46*perCent);

G4Material*PE = new G4Material("Polyethylene", 0.92*g/cm3, 2);
PE->AddElement(Hydrogen, 2);
PE->AddElement(Carbon, 4);

G4Material*shape3_mat = new G4Material("Brick", 1.65*g/cm3, 2);
shape3_mat->AddMaterial(PE, 20*perCent);
shape3_mat->AddMaterial(MMS2, 80*perCent); //p27 Sargent C6H10O2 formulas
//might add in oxygen O2 per Sargent p 27.

G4ThreeVector pos3 = G4ThreeVector(0*cm, 0*cm, 0*cm);
//sphere shape
G4double shape3_pRmin = 9*cm;
G4double shape3_pRmax = 17*cm;
G4double shape3_pSPhi = 0*deg;
G4double shape3_pDPhi = 360*deg;
G4double shape3_pSTheta = 0*deg;
G4double shape3_pDTheta = 180*deg;
G4Sphere * solidShape3 = new G4Sphere("Shape3", shape3_pRmin, shape3_pRmax, shape3_pSPhi, shape3_pDPhi, shape3_pSTheta,
shape3_pDTheta);  // size and shape half dome I hope

G4LogicalVolume* logicShape3 = 
new G4LogicalVolume(solidShape3, // its solid
			shape3_mat, // its material
			"Shape3"); //its name
new G4PVPlacement(0, pos3, logicShape3,"Shape3", logicEnv, false, 0, checkOverlaps);
fScoringVolume = logicShape3;
   
  // Shape 2
//Things that worked:
//G4double z, a, density;
//G4String name, symbol;
//G4Material* shape2_mat = new G4Material(name="Argon", z=18, a=39.95*g/mole, density=1.39*g/cm3);
//G4NistManager* Mang = G4NistManager::Instance();
//G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_CONCRETE");
//G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
//G4Material* shape2_mat = new G4Material(name="hydrogen", z=1, a=1.01*g/mole, density=1.008*g/cm3);
//Experiments


//G4Element* elO

//G4double fractionmass;

//G4Material* Na2O = nist->FindOrBuildMaterial("G4_SODIUM_MONOXIDE");
//G4Material* SiO2 = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
//G4Material* FeO = nist->FindOrBuildMaterial("G4_FERROUS_OXIDE");
//G4Material* Fe2O3 = nist->FindOrBuildMaterial("G4_FERRIC_OXIDE");
//G4Material* MgO = nist->FindOrBuildMaterial("G4_MAGNESIUM_OXIDE");
//G4Material* CaO = nist->FindOrBuildMaterial("G4_CALCIUM_OXIDE");
//G4Material* K2O = nist->FindOrBuildMaterial("G4_POTASSIUM_OXIDE");
//G4Material* TiO2 = nist->FindOrBuildMaterial("G4_TITANIUM_DIOXIDE");
//G4Element* Oxygen = nist->FindOrBuildElement("G4_O");
//G4Element* Phosphorous = nist->FindOrBuildElement("G4_P");
//G4Element* Sulfur = nist->FindOrBuildElement("G4_S");


//G4Element* Chlorine = nist->FindOrBuildElement("G4_Cl");
//G4Element* Chromium = nist->FindOrBuildElement("G4_Cr");
//G4Element* Manganese = nist->FindOrBuildElement("G4_Mn");

//G4int n;
//G4double frac;

G4Material* shape2_mat = new G4Material("Regolith", 1.7*g/cm3, 14);
shape2_mat->AddMaterial(SIO2, 46.63*perCent);
shape2_mat->AddMaterial(FeO, 12.18*perCent);
shape2_mat->AddMaterial(Fe2O3, 4.2*perCent);
shape2_mat->AddMaterial(MnO, 0.33*perCent);
shape2_mat->AddMaterial(MgO, 8.93*perCent);
shape2_mat->AddMaterial(CaO, 6.27*perCent);
shape2_mat->AddMaterial(Na2O, 3.02*perCent);
shape2_mat->AddMaterial(K2O, 0.41*perCent);
shape2_mat->AddMaterial(P2O5, 0.83*perCent);
shape2_mat->AddMaterial(Cr2O3, 0.36*perCent);
shape2_mat->AddMaterial(SO3, 4.9*perCent);
shape2_mat->AddMaterial(TiO2, 0.87*perCent);
shape2_mat->AddMaterial(Al2O3, 10.46*perCent);
shape2_mat->AddElement(Chlorine, 0.61*perCent);
//G4Material* SulfTriO = Mang->FindOrBuildMaterial("SulfurTrioxide");

//G4String name, symbol;
//G4int ncomponents;
   
  G4ThreeVector pos2 = G4ThreeVector(0*cm, -32*cm, 0*cm);
  // Trap shape       
  G4double shape2_dxa = 80.6*cm, shape2_dxb = 80.6*cm;
  G4double shape2_dya = 80.5*cm, shape2_dyb = 80.7*cm;
  G4double shape2_dz  = 10.3*cm;      
  G4Trd* solidShape2 =    
    new G4Trd("Shape2",                      //its name
              0.5*shape2_dxa, 0.5*shape2_dxb, 
              0.5*shape2_dya, 0.5*shape2_dyb, 0.5*shape2_dz); //its size
                
  G4LogicalVolume* logicShape2 =                         
   new G4LogicalVolume(solidShape2,         //its solid
                        shape2_mat,          //its material
                        "Shape2");           //its name
  G4RotationMatrix* birdo = new G4RotationMatrix();
birdo->rotateX(90*deg);
//birdo->rotateY(90*deg);
//birdo->rotateZ(25*deg);            
  new G4PVPlacement(birdo,                       //no rotation
                    pos2,                    //at position
                    logicShape2,             //its logical volume
                    "Shape2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                   checkOverlaps);          //overlaps checking


//G4Material* box_mat = nist->FindOrBuildMaterial("G4_CONCRETE");               
//G4Box* SBox = new G4Box("SogBox", 3.21*cm, 4.2*cm, 4.3*cm); //SBox size
//G4LogicalVolume* LogicSBox = new G4LogicalVolume(SBox, box_mat, "SogBox"); //logical volume
//G4RotationMatrix* Boxy = new G4RotationMatrix();
//Boxy->rotateX(39*deg);
//Boxy->rotateY(16*deg);
//Boxy->rotateZ(25*deg);
//new G4PVPlacement(Boxy, G4ThreeVector(-89,-98,0),LogicSBox,"SogBox", logicEnv, false, 0, checkOverlaps);

//
//shape 4
//
G4Material * shape4_mat = nist->FindOrBuildMaterial("G4_TISSUE_SOFT_ICRP");
G4ThreeVector pos4 = G4ThreeVector(0*cm, 0*cm, 0*cm);
G4double shape4_pRmax = 2.5*cm;
G4Orb * solidShape4 = new G4Orb("Shape4", shape4_pRmax);
G4LogicalVolume * logicShape4 = new G4LogicalVolume(solidShape4,         //its solid
                        shape4_mat,          //its material
                        "Shape4");
new G4PVPlacement(0,                       //no rotation
                    pos4,                    //at position
                    logicShape4,             //its logical volume
                    "Shape4",                //its name
                   logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                   checkOverlaps);          //overlaps checking
                
  fScoringVolume = logicShape4;

//


  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
