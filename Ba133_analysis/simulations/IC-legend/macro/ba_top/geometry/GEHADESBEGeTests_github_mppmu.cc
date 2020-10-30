//---------------------------------------------------------------------------//
//bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//                                                                           //
//                                                                           //
//                         MaGe Simulation                                   //
//                                                                           //
//      This code implementation is the intellectual property of the         //
//      MAJORANA and Gerda Collaborations. It is based on Geant4, an         //
//      intellectual property of the RD44 GEANT4 collaboration.              //
//                                                                           //
//                        *********************                              //
//                                                                           //
//    Neither the authors of this software system, nor their employing       //
//    institutes, nor the agencies providing financial support for this      //
//    work  make  any representation or  warranty, express or implied,       //
//    regarding this software system or assume any liability for its use.    //
//    By copying, distributing or modifying the Program (or any work based   //
//    on on the Program) you indicate your acceptance of this statement,     //
//    and all its terms.                                                     //
//                                                                           //
//bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//---------------------------------------------------------------------------//
//                                                          
// $Id: GEHADESBEGeTests.hh,v 1.3 2011-11-29 Bjoern $
//      
// CLASS IMPLEMENTATION:  GEHADESBEGeTests.cc
//
//---------------------------------------------------------------------------//
/**
 * SPECIAL NOTES:
 */
// 
//---------------------------------------------------------------------------//
/**
 * AUTHOR: Luciano Pandola
 * CONTACT: 
 * FIRST SUBMISSION: 
 * 
 * REVISION: 
 *  04 Jun 2009, A. di Vacri,    Revised geometry according to 
 *                                     Canberra drawings
 *  28 Apr 2010, L. Pandola,     Main revision to accomodate depGe BEGes. 
 *                               Now the basic construction via 
 *                               MGGeometryBEGeDetector is used
 *  29 Nov 2011, B. Lehnert,     Derived from GELNGSBEGeDetector class and implemented as new class in MaGe
 *  26 Jan 2012, B. Lehnert,     Adding sources
 *  20 Feb 2012, B. Lehnert,     Adding LNGS wire source
 *  09 Okt 2013, R. Falkenstein, See the following changes
 *                   03 Sep 2013 Added possibility to choose XtalMaterial via Messenger (default: EnrichedGe).
 *                               Added code for DepBEGeCryostatHolders.
 *                   07 Sep 2013 Changed material for crystal holder from pure Al to Al_alloy(EN_AW_Al_Cu6BiPb)
 *                               the support ring is unchanged.                       
 *                               Changed dimensions of HS21 encapsulation (thickness, and side-lengths) 
 *                               Changed shape of HS20 and HS29 source point, changed dimensions, changed materials.
 *                   09 Okt 2013 Changed materials for HS21 from PVC to Acrylic 
 *                               (measurement of source density of HS21 by Werner was 1.15 g/cm3, 
 *                               material is only approximation, but PVC could be excluded)
 *  10 Okt 2013, R. Falkenstein, Added plexi glass holder code from bjoern
 *                               Changed material of endcap and support ring to Al_alloy2(AL_6061_T6_AlMgSiO5)
 */
//---------------------------------------------------------------------------//
//
//#include "globals.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"
#include "G4SubtractionSolid.hh"
#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4UnionSolid.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4GeometryManager.hh"

#include "munichteststand/GEHADESBEGeTestsMessenger.hh"
#include "gerdageometry/GEGeometrySD.hh"
#include "geometry/MGGeometryGlobals.hh"
#include "geometry/MGGeometryBEGeDetector.hh"
#include "io/MGLogger.hh"

//---------------------------------------------------------------------------//

#include "munichteststand/GEHADESBEGeTests.hh"

//---------------------------------------------------------------------------//

GEHADESBEGeTests::GEHADESBEGeTests() :
		MGGeometryDetector("") {
	fMessenger = new GEHADESBEGeTestsMessenger(this);
	//Default parameters
	//Standard for now is the "CC" type detector from GELNGSBEGeDetector
	CryostatWindowThickness = 1.5 * mm;
	CryostatWallThickness = 1.5 * mm;
	CryostatDiameter = 88 * mm;
	CryostatHeight = 83 * mm;

	XtalDiameter = (2 * 37.25) * mm;
	XtalHeight = 33 * mm;
	XtalDistanceToWindow = 1.5 * mm;
	XtalDitchInnerRadius = 4.5 * mm;
	XtalDitchOuterRadius = 11.75 * mm;
	XtalDitchDepth = 2.0 * mm;
	XtalDitchOnBottom = true; // default is Ditch on bottom
	XtalCornerDiameter = 0;
	XtalCornerHeight = 0;
	XtalCornerOnBottom = false; // default is corner on top
	//raphael
	XtalMaterial = "EnrichedGe";     //default
	//end raphael
	fSpecialDetectorType = "";
	fSourceType = "";
	SourceDistance = 0;

	ActivateEnrBEGeCryostatHolders = false;
	ActivateDepBEGeCryostatHolders = false;
	ActivateHADESSourceHolder = false;

	HADESLeadCastleType = 0; //no lead castle

	fUseBRADYEnv = 0;
	fUseCollimator = 0;
	fCollimatorDistance = 1.75 * cm;
	fCollimatorPosition = 0.;

	fCryostatFillMaterial = "Vacuum"; //vaccum inside cryostat. only change via messenger for crazy RnD stuff

}

//---------------------------------------------------------------------------//

GEHADESBEGeTests::~GEHADESBEGeTests() {
	if (fMessenger)
		delete fMessenger;
}

//---------------------------------------------------------------------------//

void GEHADESBEGeTests::ConstructDetector() {

	//------------------------------------------------
	// Initialize sensitive detectors
	//------------------------------------------------
	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	G4String CrystalSDname = "/mydet/gerda/gecrystal";
	GEGeometrySD* CrystalSD = new GEGeometrySD(CrystalSDname);
	SDman->AddNewDetector(CrystalSD);

	// Retrieve materials from the Table

	G4Material* Vacuum = G4Material::GetMaterial("Vacuum");
	G4Material* Brass = G4Material::GetMaterial("Brass");
	//G4Material* ProportionalGas = G4Material::GetMaterial("ProportionalGas");
	//G4Material* N2Gas = G4Material::GetMaterial("NitrogenGas");
	G4Material* NaturalGe = G4Material::GetMaterial("NaturalGe");
	G4Material* EnrichedGe = G4Material::GetMaterial("EnrichedGe");
	G4Material* DepletedGe = G4Material::GetMaterial("DepletedGe"); //raphael
	G4Material* GeLi = G4Material::GetMaterial("Germanium/Lithium");
	G4Material* Pb = G4Material::GetMaterial("MetalLead");
	//G4Material* RSV= G4Material::GetMaterial("Steel");
	//G4Material* Fe = G4Material::GetMaterial("MetalIron");
	G4Material* Cu = G4Material::GetMaterial("MetalCopper");
	G4Material* Al = G4Material::GetMaterial("MetalAluminium");
	G4Material* Al_alloy= G4Material::GetMaterial("EN_AW_Al_Cu6BiPb"); //raphael
	G4Material* Al_alloy2= G4Material::GetMaterial("AL_6061_T6_AlMgSiO5"); //raphael
	//G4Material* Vespel = G4Material::GetMaterial("Vespel");
	G4Material* Acrylic = G4Material::GetMaterial("Acrylic");
	G4Material* PE = G4Material::GetMaterial("PE");
	G4Material* PVC = G4Material::GetMaterial("PVC");
	G4Material* Teflon = G4Material::GetMaterial("Teflon");
	G4Material* Air = G4Material::GetMaterial("Air");
	G4Material* PET = G4Material::GetMaterial("PET");
	G4Material* ssteel = G4Material::GetMaterial("Steel");       // moje
	G4Material* water = G4Material::GetMaterial("Water");       // moje
	G4Material* PS = G4Material::GetMaterial("PS");       //raphael (Polystyrene)

	G4Material* AmericiumOxide = G4Material::GetMaterial("AmericiumOxide");

	//G4Material* Quartz = G4Material::GetMaterial("Quartz");

	G4VisAttributes * Gecolour = new G4VisAttributes(
			G4Colour(255 / 255., 1 / 255., 1 / 255.));  //red
	G4VisAttributes * Cucolour = new G4VisAttributes(
			G4Colour(1 / 255., 255 / 255., 33 / 255.));  //green !
	G4VisAttributes * Alcolour = new G4VisAttributes(
			G4Colour(1 / 255., 1 / 255., 255 / 255.));   //blue
	//raphael
	G4VisAttributes * Alalloycolour = new G4VisAttributes(
			G4Colour(193 / 255., 193 / 255., 255 / 255.));   //pale blue
	//end raphael
	G4VisAttributes * PEcolour = new G4VisAttributes(
			G4Colour(226 / 255., 163 / 255., 29 / 255.));  //orange
	G4VisAttributes * Vacuumcolor = new G4VisAttributes(
			G4Colour(210 / 255., 210 / 255., 210 / 255.)); //light blue-grey , moje

	//G4VisAttributes * ProportionalGascolour= new G4VisAttributes( G4Colour(154/255. ,237/255. ,193/255. ));  //blue-green
	G4VisAttributes * Pbcolour = new G4VisAttributes(
			G4Colour(171 / 255., 171 / 255., 195 / 255.));  //gray
	G4VisAttributes * GeLicolour = new G4VisAttributes(
			G4Colour(187 / 255., 28 / 255., 0 / 255.)); //almost black ( -> 'dead'...)
	//G4VisAttributes * N2Gascolour = new G4VisAttributes( G4Colour(150/255. ,150/255. ,255/255. ));   //light blue
	G4VisAttributes * Tefloncolour = new G4VisAttributes(
			G4Colour(150 / 255., 213 / 255., 150 / 255.));
	G4VisAttributes * PETcolour = new G4VisAttributes(
			G4Colour(255 / 255., 71 / 255., 33 / 255.));  //red

	G4VisAttributes * Acryliccolour = new G4VisAttributes(
			G4Colour(231 / 255., 217 / 255., 240 / 255.));  //light blue-grey
	//G4VisAttributes * RSVcolour= new G4VisAttributes( G4Colour(50/255. ,150/255. ,50/255. ));
	G4VisAttributes * Brasscolour = new G4VisAttributes(
			G4Colour(226 / 255., 163 / 255., 29 / 255.));  //orange-brown

	G4VisAttributes * Ceramiccolour = new G4VisAttributes(
			G4Colour(255 / 255., 71 / 255., 33 / 255.));  //red , moje
	G4VisAttributes * Ssteelcolour = new G4VisAttributes(
			G4Colour(231 / 255., 217 / 255., 240 / 255.)); //light blue-grey , moje

	G4VisAttributes * AmOxcolour = new G4VisAttributes(
			G4Colour(209 / 255., 134 / 255., 0 / 255.));  //orange

	// ListMaterials
	// G4NistManager* manager = G4NistManager::Instance();
	// manager->ListMaterials("all");

	// ####################################################################
	// ------------------- experimental hall (world volume) ---------------

	G4double expHall_x = 2 * m;
	G4double expHall_y = 2 * m;
	G4double expHall_z = 2 * m;
	G4Box* experimentalHall = new G4Box("experimentalHall", expHall_x,
			expHall_y, expHall_z);
	G4LogicalVolume* experimentalHallLogical = new G4LogicalVolume(
			experimentalHall, Air, "experimentalHallLogical", 0, 0, 0);
	experimentalHallLogical->SetVisAttributes(G4VisAttributes::Invisible);

	// ####################################################################
	// ------------------- detector housing -------------------------------

	//Cryostat shell (made of Al)
	//raphael: changed material from pure aluminum to Al_alloy2
	G4Tubs* cryostatWall = new G4Tubs("cryostatWall", 0, CryostatDiameter / 2,
			CryostatHeight / 2., 0, 360. * deg);    //Al solid cylinder
	G4LogicalVolume *cryostatWallLogical = new G4LogicalVolume(cryostatWall, Al_alloy2,
			"cryostatWallLogical");
	cryostatWallLogical->SetVisAttributes(Alcolour);
	new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), cryostatWallLogical,
			"cryostatWallPhysical", experimentalHallLogical, false, 0);

	//Vacuum inside the Cryostat shell
	G4Tubs* cryostatVacuum = new G4Tubs("cryostatVacuum", 0.,
			CryostatDiameter / 2 - CryostatWallThickness,
			CryostatHeight / 2. - CryostatWindowThickness, 0, 360. * deg);
	G4LogicalVolume* cryostatVacuumLogical = new G4LogicalVolume(cryostatVacuum,
			G4Material::GetMaterial(fCryostatFillMaterial), "cryostatVacuumLogical", 0, 0, 0);
	cryostatVacuumLogical->SetVisAttributes(G4VisAttributes::Invisible);
	new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), cryostatVacuumLogical,
			"cryostatVacuumPhysical", cryostatWallLogical, false, 0);

	cryostatVacuumLogical->SetVisAttributes(Vacuumcolor);

	// ####################################################################
	// ------------------- detector without dead layer --------------------

	if (fSpecialDetectorType == "") //no special geometry, construct BEGe from template class
			{
		MGGeometryBEGeDetector* BEGeTemplate = new MGGeometryBEGeDetector(
				"logicDet");

		BEGeTemplate->SetRadius(XtalDiameter / 2.);
		BEGeTemplate->SetHeight(XtalHeight);
		BEGeTemplate->SetDitchInnerRadius(XtalDitchInnerRadius);
		BEGeTemplate->SetDitchOuterRadius(XtalDitchOuterRadius);
		BEGeTemplate->SetDitchDepth(XtalDitchDepth);
		//raphael
		BEGeTemplate->SetG4MaterialName(XtalMaterial);
		//end raphael
		BEGeTemplate->SetDitchBelow(XtalDitchOnBottom);
		BEGeTemplate->SetCornerDiameter(XtalCornerDiameter);
		BEGeTemplate->SetCornerHeight(XtalCornerHeight);
		BEGeTemplate->SetCornerOnBottom(XtalCornerOnBottom);

		BEGeTemplate->ConstructDetector();
		G4LogicalVolume *BEGeLogical = BEGeTemplate->GetDetectorLogical();
		BEGeLogical->SetVisAttributes(Gecolour);

		// calculating BEGe center z position with cryostat hight and XtalDistanceToWindow
		G4double zPositionBEGeCenter = CryostatHeight / 2.
				- CryostatWindowThickness - XtalDistanceToWindow
				- XtalHeight / 2.;

		new G4PVPlacement(0, G4ThreeVector(0, 0, zPositionBEGeCenter),
				BEGeLogical, "BEGeTemplatePhysical", cryostatVacuumLogical,
				false, 0);

		BEGeLogical->SetSensitiveDetector(CrystalSD);

	} else if (fSpecialDetectorType == "some geometry") //special geometry
			{

		// implement special geometries
		MGLog(error)<<" This geometry is not yet implemented" << endlog;

	}
	else
	{
		MGLog(error) <<fSpecialDetectorType <<" geometry not implemented" << endlog;
	}

	SetDetectorLogical(experimentalHallLogical);
	SetDetectorName("ExperimentalHall");

	// ####################################################################
	// ------------------- holders ----------------------------------------
	G4double zPositionBEGeCenter = CryostatHeight / 2 - CryostatWindowThickness
			- XtalDistanceToWindow - XtalHeight / 2;

	if (ActivateEnrBEGeCryostatHolders == true) {

		//HD1000 part (xtal cup):
		// approximated with high density PE (94g/cm3)

		G4double cupSeg1zPos = zPositionBEGeCenter + XtalHeight / 2 + 1.5 * mm
				- 55. * mm / 2;
		G4Tubs* cupSeg1 = new G4Tubs("cupSeg1",
				XtalDiameter / 2 + 0.5 * mm / 2.,
				XtalDiameter / 2. + 1.5 * mm / 2, 55. * mm / 2, 0, 360. * deg);
		G4LogicalVolume *cupSeg1Logical = new G4LogicalVolume(cupSeg1, PE,
				"cupSeg1Logical");
		cupSeg1Logical->SetVisAttributes(PEcolour);
		new G4PVPlacement(0, G4ThreeVector(0., 0., cupSeg1zPos), cupSeg1Logical,
				"cupSeg1Physical", cryostatVacuumLogical, false, 0);

		G4double cupSeg2zPos = zPositionBEGeCenter + XtalHeight / 2 + 0.5 * mm
				+ 0.5 * mm;
		G4Tubs* cupSeg2 = new G4Tubs("cupSeg2", 0. * mm / 2.,
				XtalDiameter / 2. + 1.5 * mm / 2, 1. * mm / 2, 0, 360. * deg);
		G4LogicalVolume *cupSeg2Logical = new G4LogicalVolume(cupSeg2, PE,
				"cupSeg2Logical");
		cupSeg2Logical->SetVisAttributes(PEcolour);
		new G4PVPlacement(0, G4ThreeVector(0., 0., cupSeg2zPos), cupSeg2Logical,
				"cupSeg2Physical", cryostatVacuumLogical, false, 0);

		//Aluminum part (Holder)
		//raphael: changed material from pure aluminum to aluminum alloy
		G4double holderSeg1zPos = zPositionBEGeCenter + XtalHeight / 2
				+ 1.5 * mm - 62.5 * mm / 2; // 1mm due to xtal cub
		G4Tubs* holderSeg1 = new G4Tubs("holderSeg1",
				XtalDiameter / 2. + 2. * mm / 2,
				XtalDiameter / 2. + 5. * mm / 2, 62.5 * mm / 2, 0, 360. * deg);
		G4LogicalVolume *holderSeg1Logical = new G4LogicalVolume(holderSeg1, Al_alloy,
				"holderSeg1Logical"); //raphael
		holderSeg1Logical->SetVisAttributes(Alalloycolour); //raphael
		new G4PVPlacement(0, G4ThreeVector(0., 0., holderSeg1zPos),
				holderSeg1Logical, "holderSeg1Physical", cryostatVacuumLogical,
				false, 0);

		G4double holderSeg2zPos = holderSeg1zPos - 62.5 * mm / 2 - 3. * mm / 2;
		G4Tubs* holderSeg2 = new G4Tubs("holderSeg2", 53. * mm / 2,
				XtalDiameter / 2. + 5. * mm / 2, 3. * mm / 2, 0, 360. * deg);
		G4LogicalVolume *holderSeg2Logical = new G4LogicalVolume(holderSeg2, Al_alloy,
				"holderSeg2Logical"); //raphael
		holderSeg2Logical->SetVisAttributes(Alalloycolour); //raphael
		new G4PVPlacement(0, G4ThreeVector(0., 0., holderSeg2zPos),
				holderSeg2Logical, "holderSeg2Physical", cryostatVacuumLogical,
				false, 0);

		G4double holderSeg3zPos = holderSeg2zPos - 3. * mm / 2 - 12. * mm / 2;
		G4Tubs* holderSeg3 = new G4Tubs("holderSeg3", 53. * mm / 2,
				60. * mm / 2, 12. * mm / 2, 0, 360. * deg);
		G4LogicalVolume *holderSeg3Logical = new G4LogicalVolume(holderSeg3, Al_alloy,
				"holderSeg3Logical"); //raphael
		holderSeg3Logical->SetVisAttributes(Alalloycolour); //raphael
		new G4PVPlacement(0, G4ThreeVector(0., 0., holderSeg3zPos),
				holderSeg3Logical, "holderSeg3Physical", cryostatVacuumLogical,
				false, 0);
		//raphael: changed material from pure aluminum to Al_alloy2
		G4double holderSeg4zPos = holderSeg1zPos + 62.5 * mm / 2 - 19. * mm;
		G4Tubs* holderSeg4 = new G4Tubs("holderSeg4",
				XtalDiameter / 2. + 5. * mm / 2,
				XtalDiameter / 2. + 8. * mm / 2, 8.6 * mm / 2, 0, 360. * deg);
		G4LogicalVolume *holderSeg4Logical = new G4LogicalVolume(holderSeg4, Al_alloy2,
				"holderSeg4Logical");
		holderSeg4Logical->SetVisAttributes(Alcolour);
		new G4PVPlacement(0, G4ThreeVector(0., 0., holderSeg4zPos),
				holderSeg4Logical, "holderSeg4Physical", cryostatVacuumLogical,
				false, 0);

		//Copper part (Electronic base):

		G4double baseSeg1zPos = holderSeg1zPos - 62.5 * mm / 2 - 17. * mm
				- 26. * mm / 2;
		G4Tubs* baseSeg1 = new G4Tubs("baseSeg1", 0. * mm / 2, 26. * mm / 2,
				26. * mm / 2, 0, 360. * deg);

		// cut out upper interior tube of baseSeg1
		G4Tubs* baseNegSeg1 = new G4Tubs("baseNegSeg1", 0. * mm / 2,
				8. * mm / 2, 13. * mm / 2, 0, 360. * deg);
		G4VSolid* baseSeg1_1 = new G4SubtractionSolid("baseSeg1_1", baseSeg1,
				baseNegSeg1, 0, G4ThreeVector(0, 0, +13 / 2.));

		// cut out lower interior tube of baseSeg1
		G4Tubs* baseNegSeg2 = new G4Tubs("baseNegSeg2", 0. * mm / 2,
				3. * mm / 2, 15. * mm / 2, 0, 360. * deg);
		G4VSolid* baseSeg1_2 = new G4SubtractionSolid("baseSeg1_2", baseSeg1_1,
				baseNegSeg2, 0, G4ThreeVector(0, 0, -13 / 2.));

		// cut out exterior of basSeg1_1
		G4Tubs* baseNegSeg3 = new G4Tubs("baseNegSeg3", 19. * mm / 2,
				27. * mm / 2, 11. * mm / 2, 0, 360. * deg);
		G4VSolid* baseSeg1_3 = new G4SubtractionSolid("baseSeg1_3", baseSeg1_2,
				baseNegSeg3, 0, G4ThreeVector(0, 0, +11 / 2. + 2));

		G4LogicalVolume *baseSeg1Logical = new G4LogicalVolume(baseSeg1_3, Cu,
				"baseSeg1Logical");
		baseSeg1Logical->SetVisAttributes(Cucolour);

		new G4PVPlacement(0, G4ThreeVector(0., 0., baseSeg1zPos),
				baseSeg1Logical, "baseSeg1Physical", cryostatVacuumLogical,
				false, 0);

		G4double baseSeg2zPos = holderSeg1zPos - 62.5 * mm / 2 - 17. * mm
				+ 3. * mm / 2;
		G4Tubs* baseSeg2 = new G4Tubs("baseSeg2", 8. * mm / 2, 38. * mm / 2,
				3. * mm / 2, 0, 360. * deg);
		G4LogicalVolume *baseSeg2Logical = new G4LogicalVolume(baseSeg2, Cu,
				"baseSeg2Logical");
		baseSeg2Logical->SetVisAttributes(Cucolour);
		new G4PVPlacement(0, G4ThreeVector(0., 0., baseSeg2zPos),
				baseSeg2Logical, "baseSeg2Physical", cryostatVacuumLogical,
				false, 0);

		G4double baseSeg3zPos = holderSeg1zPos - 62.5 * mm / 2 - 17. * mm / 2;
		G4Tubs* baseSeg3 = new G4Tubs("baseSeg3", 38. * mm / 2, 52.8 * mm / 2,
				17. * mm / 2, 0, 360. * deg);
		G4LogicalVolume *baseSeg3Logical = new G4LogicalVolume(baseSeg3, Cu,
				"baseSeg3Logical");
		baseSeg3Logical->SetVisAttributes(Cucolour);
		new G4PVPlacement(0, G4ThreeVector(0., 0., baseSeg3zPos),
				baseSeg3Logical, "baseSeg3Physical", cryostatVacuumLogical,
				false, 0);

		G4double baseSeg4zPos = baseSeg2zPos - 3. * mm / 2 + 2. * mm / 2;
		G4Tubs* baseSeg4 = new G4Tubs("baseSeg4", 52.8 * mm / 2, 60. * mm / 2,
				2. * mm / 2, 0, 360. * deg);
		G4LogicalVolume *baseSeg4Logical = new G4LogicalVolume(baseSeg4, Cu,
				"baseSeg4Logical");
		baseSeg4Logical->SetVisAttributes(Cucolour);
		new G4PVPlacement(0, G4ThreeVector(0., 0., baseSeg4zPos),
				baseSeg4Logical, "baseSeg4Physical", cryostatVacuumLogical,
				false, 0);

	} // end if(ActivateEnrBEGeCryostatHolders=true)
	else if(ActivateDepBEGeCryostatHolders == true)
	{
	  //raphael (with information from Katharina)
	  // TUBEGE COPPER HOLDER

	    // CONSTANTS DEFINED ALL POSITIVE!
	    G4double mantelHeight = 52.*mm;
	    G4double mantelThickness = 1.5*mm;
	    G4double distXtalMantel = 1.*mm;
	    G4double mantelInnerDiameter = XtalDiameter + 2.*distXtalMantel;
	    G4double mantelOuterDiameter = mantelInnerDiameter + 2.*mantelThickness;
	    G4double mantelShiftToXtalTop = 1.5*mm;
	    G4double firstRingShiftToMantelTop = 1.5*mm;
	    G4double firstRingHeight = 7.*mm;
	    G4double firstRingThickness = 1.5*mm;
	    G4double secondRingShiftToMantelTop = 13.*mm;
	    G4double secondRingHeight = firstRingHeight;
	    G4double secondRingThickness = firstRingThickness;
	    G4double holderBaseHeight = 10.*mm;
	    G4double holderBaseDiameter = mantelOuterDiameter;
	    
	    // MANTEL
	    G4double holderMantelPos = zPositionBEGeCenter + XtalHeight/2. - mantelShiftToXtalTop - mantelHeight/2.;
	    G4Tubs* holderMantel
		= new G4Tubs("holderMantel", mantelInnerDiameter/2., mantelOuterDiameter/2., mantelHeight/2., 0, 360.*deg);
	    G4LogicalVolume *holderMantelLogical
		= new G4LogicalVolume(holderMantel, Cu, "holderMantelLogical");  // CREATE LOGICAL VOLUME
	    holderMantelLogical->SetVisAttributes(Cucolour);
	    new G4PVPlacement(0, G4ThreeVector(0.,0.,holderMantelPos), holderMantelLogical, "holderMantelPhysical", cryostatVacuumLogical, false, 0);
	    
	    // FIRST RING
	    G4double holderRing1Pos = holderMantelPos + mantelHeight/2. - firstRingShiftToMantelTop - firstRingHeight/2.;
	    G4Tubs* holderRing1
		= new G4Tubs("holderRing1", mantelOuterDiameter/2., mantelOuterDiameter/2. + firstRingThickness, firstRingHeight/2., 0, 360.*deg);
	    G4LogicalVolume *holderRing1Logical
		= new G4LogicalVolume(holderRing1, Cu, "holderRing1Logical");  // CREATE LOGICAL VOLUME
	    holderRing1Logical->SetVisAttributes(Cucolour);
	    new G4PVPlacement(0, G4ThreeVector(0.,0.,holderRing1Pos), holderRing1Logical, "holderRing1Physical", cryostatVacuumLogical, false, 0);
	    
	    // SECOND RING
	    G4double holderRing2Pos = holderMantelPos + mantelHeight/2. - secondRingShiftToMantelTop - secondRingHeight/2.;
	    G4Tubs* holderRing2
		= new G4Tubs("holderRing2", mantelOuterDiameter/2., mantelOuterDiameter/2. + secondRingThickness, firstRingHeight/2., 0, 360.*deg);
	    G4LogicalVolume *holderRing2Logical
		= new G4LogicalVolume(holderRing2, Cu, "holderRing2Logical");  // CREATE LOGICAL VOLUME
	    holderRing2Logical->SetVisAttributes(Cucolour);
	    new G4PVPlacement(0, G4ThreeVector(0.,0.,holderRing2Pos), holderRing2Logical, "holderRing2Physical", cryostatVacuumLogical, false, 0);
	    
	    // BASE
	    G4double holderBasePos = holderMantelPos - mantelHeight/2. - holderBaseHeight/2.;
	    G4Tubs* holderBase
		= new G4Tubs("holderBase", 0, holderBaseDiameter/2., holderBaseHeight/2., 0, 360.*deg);
	    G4LogicalVolume *holderBaseLogical
		= new G4LogicalVolume(holderBase, Cu, "holderBaseLogical");  // CREATE LOGICAL VOLUME
	    holderBaseLogical->SetVisAttributes(Cucolour);
	    new G4PVPlacement(0, G4ThreeVector(0.,0.,holderBasePos), holderBaseLogical, "holderBasePhysical", cryostatVacuumLogical, false, 0); 

	    //PE-cup 
	    //still missing
	    
	    //end raphael
	}

	// by bjoern
	// ##############################################################################
	// construct acrylic source holder composed of independent generic volumes parts
	// construct only if messenger is set true (default=false)
	// top lid of source holder will be constructed with the source to avoid overlaps
	// the top lid is only implemented for HS21 and PTBPointlike
	
	
	G4double SHThickness = 21./2.*mm; // from technical drawing
	G4double SHOuterRadius = 108./2.*mm; // from technical drawing
	G4double SHInnerRadius = SHOuterRadius - SHThickness;
	
	G4double SHLength = SourceDistance+3*mm; // from source distance and source thickness
	G4double SHBottomExtensionHeight = 4.6*mm;    //raphael: according to Katharinas and Raphaels measurement it should be 4.7
	G4double SHLidHeight = 8*mm;
	
	
	G4ThreeVector* SHVector = new G4ThreeVector(0,0,CryostatHeight/2. + SHLength/2.);
	
	if (ActivateHADESSourceHolder == true) {
	  
	  
	  
	  // main SH volume
	  G4Tubs* SHWall = new G4Tubs("SHWall", SHInnerRadius, SHOuterRadius,
				      SHLength / 2., 0, 360. * deg);
	  G4LogicalVolume *SHWallLogical = new G4LogicalVolume(SHWall, Acrylic,
							       "SHWallLogical");
	  SHWallLogical->SetVisAttributes(Acryliccolour);
	  new G4PVPlacement(0, *SHVector, SHWallLogical,
			    "SHWallPhysical", experimentalHallLogical, false, 0);
	  
	  
	  // small SH extension around endcap
	  G4Tubs* SHBottomExtension = new G4Tubs("SHBottomExtension", 102./2., SHOuterRadius,
						 SHBottomExtensionHeight / 2., 0, 360. * deg);
	  G4LogicalVolume *SHBottomExtensionLogical = new G4LogicalVolume(SHBottomExtension, Acrylic,
									  "SHBottomExtensionLogical");
	  SHBottomExtensionLogical->SetVisAttributes(Acryliccolour);
	  new G4PVPlacement(0, *SHVector+G4ThreeVector(0,0,-SHLength/2.-SHBottomExtensionHeight/2.), SHBottomExtensionLogical,
			    "SHBottomExtensionPhysical", experimentalHallLogical, false, 0);
	}
	//end by bjoern

	
	// ####################################################################
	// ------------------- sources ----------------------------------------

	//#################################################
	//Tueb source: Information by Katharinas
	if (fSourceType == "Tueb") {
		MGLog(trace)<<"Loading "<<fSourceType<<" source... "<<endlog;

		// Validation Measurement in Tuebingen

		//     G4double SourceCoatingDistanceToWindow = SourceDistance;
		G4double SourceCoatingHeight = 2*mm;
		G4double SourceCoatingLength = 2*cm;
		G4double SourceCoatingWidth = 1*cm;

		G4double SourcePointRadius = .05*mm;

		G4ThreeVector* SourceVector = new G4ThreeVector(0,0,CryostatHeight/2 + SourceCoatingHeight/2 + SourceDistance);

		//source coating
		G4Box *SourceCoating
		= new G4Box("SourceCoating",SourceCoatingLength/2., SourceCoatingWidth/2., SourceCoatingHeight/2.);
		G4LogicalVolume* SourceCoatingLogical
		= new G4LogicalVolume(SourceCoating,Acrylic , "SourceCoatingLogical");
		SourceCoatingLogical->SetVisAttributes(Acryliccolour);

		new G4PVPlacement(0, *SourceVector, SourceCoatingLogical,"SourceCoatingPhysical",experimentalHallLogical, false, 0);

		//source volume
		G4Sphere *Source
		= new G4Sphere("Source",0, SourcePointRadius,0, 360.*deg,0, 180.*deg);
		G4LogicalVolume* SourceLogical
		= new G4LogicalVolume(Source,Acrylic , "SourceLogical");
		SourceLogical->SetVisAttributes(PETcolour);

		new G4PVPlacement(0, G4ThreeVector(0,0,0), SourceLogical,"SourcePhysical",SourceCoatingLogical, false, 0);

	}

	//#################################################
	//PTB source: Information by Erica
	else if(fSourceType=="PTBPointlike")
	{
		MGLog(trace)<<"Loading "<<fSourceType<<" source... "<<endlog;

		/*
		 source thickness according to specs:
		 mass per area: 21.3+-1.8 mg*cm-2
		 density PE (wiki): 0.915-0.97 g/cm3
		 thickness: 0.201-0.253 mm
		 average thickness: 0.227mm

		 The source is constructed with an outer Al ring that holds two layers of thin PE foil that are sandwitching a "very thin compact-grain 5mm layer of source material"

		 */
		G4double fSourceThickness = 0.227 *mm; // thickness of one side of foil
		G4double fSourceRingThickness = 3. *mm;// thickness of one side of foil

		//source location
		G4ThreeVector* SourceVector = new G4ThreeVector(0,0,CryostatHeight/2. + fSourceRingThickness/2. + SourceDistance);

		G4RotationMatrix* SourceRotation = new G4RotationMatrix();

		//source ring
		G4Tubs* SourceRing = new G4Tubs("SourceRing",
				20.001/2. * mm,
				30./2. * mm,
				fSourceRingThickness/2.,
				0.0 * deg,
				360 * deg);
		G4LogicalVolume* SourceRingLogical
		= new G4LogicalVolume(SourceRing, Al, "SourceRingLogical");
		new G4PVPlacement(SourceRotation, *SourceVector, SourceRingLogical,"SourceRingPhysical",experimentalHallLogical, false, 0);

		// source foil
		G4Tubs* SourceFoil = new G4Tubs("SourceFoil",
				0. * mm,
				20./2. * mm,
				fSourceThickness,
				0.0 * deg,
				360 * deg);
		G4LogicalVolume* SourceFoilLogical
		= new G4LogicalVolume(SourceFoil,PE , "SourceFoilLogical");
		new G4PVPlacement(0, *SourceVector , SourceFoilLogical,"SourceFoilPhysical",experimentalHallLogical, false, 0);

		// source point
		G4Tubs* SourcePoint = new G4Tubs("SourcePoint",
				0. * mm,
				5./2. * mm,
				0.005 * mm,
				0.0 * deg,
				360 * deg);
		G4LogicalVolume* SourcePointLogical
		= new G4LogicalVolume(SourcePoint, PE , "SourcePointLogical");
		SourcePointLogical->SetVisAttributes(PETcolour);
		new G4PVPlacement(0, G4ThreeVector(0,0,0) , SourcePointLogical,"SourcePhysical",SourceFoilLogical, false, 0);

		// by bjoern
		// constructing source holder top plate
		if (ActivateHADESSourceHolder == true) {

		  G4double SHLidInnerRadius = 15.0*mm + 0.1*mm; // radius of outer ring of source + epsilon
		  
		  G4Tubs* SHLid = new G4Tubs("SHLid", SHLidInnerRadius, SHInnerRadius,
					     SHLidHeight / 2., 0, 360. * deg);
		  G4LogicalVolume *SHLidLogical = new G4LogicalVolume(SHLid, Acrylic,
								      "SHLidLogical");
		  SHLidLogical->SetVisAttributes(Acryliccolour);
		  new G4PVPlacement(0, *SourceVector + G4ThreeVector(0,0,fSourceRingThickness/2.-SHLidHeight/2.), SHLidLogical,
				    "SHLidPhysical", experimentalHallLogical, false, 0);
		}
		//end by bjoern
	}
	
	//#################################################
	//AlPill source: Information by Erica
	else if((fSourceType=="AlPillUR")||(fSourceType=="AlPillUD")||(fSourceType=="AlPill"))
	{
		// asymmetric pill with 1.2mm Al on top, 1.5mm on bottom and 0.3mm cavity

		MGLog(trace)<<"Loading "<<fSourceType<<" source... "<<endlog;

		G4double fSourceThickness = 3.*mm;// thickness of the pill

		//source location
		G4ThreeVector* SourceVector = new G4ThreeVector(0,0,CryostatHeight/2. + fSourceThickness/2. + SourceDistance);


		G4RotationMatrix* SourceRotation = new G4RotationMatrix();

		// put upside down if AlPillUD else, leave itÊ
		if (fSourceType=="AlPillUD")
			{
				SourceRotation ->rotateY(180*deg);
			}

		//source pill
		G4Tubs* SourcePill = new G4Tubs("SourcePill",
				0. * mm,
				25./2. * mm,
				fSourceThickness/2.,
				0.0 * deg,
				360 * deg);
		G4LogicalVolume* SourcePillLogical
		= new G4LogicalVolume(SourcePill, Al, "SourcePillLogical");
		new G4PVPlacement(SourceRotation, *SourceVector, SourcePillLogical,"SourcePillPhysical",experimentalHallLogical, false, 0);

		// source point
		G4Tubs* SourcePoint = new G4Tubs("SourcePoint",
				0. * mm,
				15./2. * mm,
				0.3/2 * mm,
				0.0 * deg,
				360 * deg);
		G4LogicalVolume* SourcePointLogical
		= new G4LogicalVolume(SourcePoint, water , "SourcePointLogical");
		SourcePointLogical->SetVisAttributes(PETcolour);
		new G4PVPlacement(0, G4ThreeVector(0,0,+0.15*mm) , SourcePointLogical,"SourcePhysical",SourcePillLogical, false, 0);
	}

	//#################################################
	//HS29 Ba133 source: Infomation by Werner

	else if((fSourceType=="HS29")||(fSourceType=="HS29like"))
	{
		// symmetric rectangular shape 2.26mm thickness with 1mm diameter source  //raphael (changed from 2.28 to 2.26)

		MGLog(trace)<<"Loading "<<fSourceType<<" source... "<<endlog;

		G4double fSourceThickness = 2.26*mm;// thickness of the pill //raphael (changed from 2.28 to 2.26)

		//source location
		G4ThreeVector* SourceVector = new G4ThreeVector(0,0,CryostatHeight/2. + fSourceThickness/2. + SourceDistance);

		G4RotationMatrix* SourceRotation = new G4RotationMatrix();

		//source frame
		G4Box *SourceFrame = new G4Box("SourceFrame",
				10./2. * mm,
				23./2. * mm,
				fSourceThickness/2.);

		G4LogicalVolume* SourceFrameLogical
		= new G4LogicalVolume(SourceFrame, PS, "SourceFrameLogical"); //raphael (changed from PVC to PS)
		new G4PVPlacement(SourceRotation, *SourceVector, SourceFrameLogical,"SourceFramePhysical",experimentalHallLogical, false, 0);

		// source point
		//raphael: change from cylinder (2mm dia/1mm height), to sphere (1mm dia) 
		
		//G4Tubs* SourcePoint = new G4Tubs("SourcePoint",
		//		0. * mm,
		//		2./2. * mm,
		//		1./2 * mm,
		//		0.0 * deg,
		//		360 * deg);
		
		G4Sphere* SourcePoint = new G4Sphere("SourcePoint",
				0. * mm,
				1./2. * mm,
				0.0 * deg,
				360 * deg,
				0.0 * deg,
				180 * deg);
		//end raphael
		G4LogicalVolume* SourcePointLogical
		= new G4LogicalVolume(SourcePoint, PS, "SourcePointLogical"); //raphael (changed from PVC to PS)
		SourcePointLogical->SetVisAttributes(PETcolour);
		new G4PVPlacement(0, G4ThreeVector(0,0,0) , SourcePointLogical,"SourcePhysical",SourceFrameLogical, false, 0);

		// raphael
		// constructing source holder top plate
		if (ActivateHADESSourceHolder == true) {
		  
		  G4double SHLidInnerRadius = pow(pow(10.,2)+pow(23.,2),0.5)/2.*mm + 0.1*mm; // diagonal of source frame rectangular + epsilon 

		  G4Tubs* SHLid = new G4Tubs("SHLid", SHLidInnerRadius, SHInnerRadius,
					     SHLidHeight / 2., 0, 360. * deg);
		  G4LogicalVolume *SHLidLogical = new G4LogicalVolume(SHLid, Acrylic,
								      "SHLidLogical");
		  SHLidLogical->SetVisAttributes(Acryliccolour);
		  new G4PVPlacement(0, *SourceVector + G4ThreeVector(0,0,fSourceThickness/2.-SHLidHeight/2.), SHLidLogical,
				    "SHLidPhysical", experimentalHallLogical, false, 0);
		}
		// end raphael
	}

	//#################################################
	//HS20 Am241 source: Infomation by Werner
	else if((fSourceType=="HS20")||(fSourceType=="HS20like"))
	{
		// symmetric rectangular shape 2.32mm thickness with 1mm diameter source //raphael changed from 2.3 to 2.32

		MGLog(trace)<<"Loading "<<fSourceType<<" source... "<<endlog;

		G4double fSourceThickness = 2.32*mm;// thickness of the pill  //raphael changed from 2.3 to 2.32

		//source location
		G4ThreeVector* SourceVector = new G4ThreeVector(0,0,CryostatHeight/2. + fSourceThickness/2. + SourceDistance);

		G4RotationMatrix* SourceRotation = new G4RotationMatrix();

		//source frame
		G4Box *SourceFrame = new G4Box("SourceFrame",
				10./2. * mm,
				23./2. * mm,
				fSourceThickness/2.);

		G4LogicalVolume* SourceFrameLogical
		  = new G4LogicalVolume(SourceFrame, PS, "SourceFrameLogical"); //raphael (changed from PVC to PS)
		new G4PVPlacement(SourceRotation, *SourceVector, SourceFrameLogical,"SourceFramePhysical",experimentalHallLogical, false, 0);

		// source point
		//raphael: change from cylinder (2mm dia/1mm height), to sphere (1mm dia) 
		/*
		G4Tubs* SourcePoint = new G4Tubs("SourcePoint",
				0. * mm,
				2./2. * mm,
				1./2 * mm,
				0.0 * deg,
				360 * deg);
		*/
		G4Sphere* SourcePoint = new G4Sphere("SourcePoint",
				0. * mm,
				1./2. * mm,
				0.0 * deg,
				360 * deg,
				0.0 * deg,
				180 * deg);
		//end raphael
		G4LogicalVolume* SourcePointLogical
		= new G4LogicalVolume(SourcePoint, PS , "SourcePointLogical");  //raphael (changed from PVC to PS)
		SourcePointLogical->SetVisAttributes(PETcolour);
		new G4PVPlacement(0, G4ThreeVector(0,0,0) , SourcePointLogical,"SourcePhysical",SourceFrameLogical, false, 0);

		// raphael
		// constructing source holder top plate
		if (ActivateHADESSourceHolder == true) {
		  
		  G4double SHLidInnerRadius = pow(pow(10.,2)+pow(23.,2),0.5)/2.*mm + 0.1*mm; // diagonal of source frame rectangular + epsilon
		  
		  G4Tubs* SHLid = new G4Tubs("SHLid", SHLidInnerRadius, SHInnerRadius,
					     SHLidHeight / 2., 0, 360. * deg);
		  G4LogicalVolume *SHLidLogical = new G4LogicalVolume(SHLid, Acrylic,
								      "SHLidLogical");
		  SHLidLogical->SetVisAttributes(Acryliccolour);
		  new G4PVPlacement(0, *SourceVector + G4ThreeVector(0,0,fSourceThickness/2.-SHLidHeight/2.), SHLidLogical,
				    "SHLidPhysical", experimentalHallLogical, false, 0);
		}
		//end raphael
	}

	//#################################################
	//HS21 Am241 source: Infomation meassured
	else if((fSourceType=="HS21")||(fSourceType=="HS21like"))
	{
		// symmetric rectangular shape 2.02mm thickness with 1mm spherical source //raphael (changed from 2.30 to 2.02)

		MGLog(trace)<<"Loading "<<fSourceType<<" source... "<<endlog;

		G4double fSourceThickness = 2.02*mm; //raphael (changed from 2.30 to 2.02)

		//source location
		G4ThreeVector* SourceVector = new G4ThreeVector(0,0,CryostatHeight/2. + fSourceThickness/2. + SourceDistance);

		G4RotationMatrix* SourceRotation = new G4RotationMatrix();

		//source frame
		G4Box *SourceFrame = new G4Box("SourceFrame",
				11.08/2. * mm,
				23.48/2. * mm,
				fSourceThickness/2.); //raphael (changed dimensions)

		G4LogicalVolume* SourceFrameLogical
		= new G4LogicalVolume(SourceFrame, Acrylic, "SourceFrameLogical");  //raphael (changed material from PVC to Acrylic)
		new G4PVPlacement(SourceRotation, *SourceVector, SourceFrameLogical,"SourceFramePhysical",experimentalHallLogical, false, 0);

		// source point
		G4Sphere* SourcePoint = new G4Sphere("SourcePoint",
				0. * mm,
				1./2. * mm,
				0.0 * deg,
				360 * deg,
				0.0 * deg,
				180 * deg);

		G4LogicalVolume* SourcePointLogical
		= new G4LogicalVolume(SourcePoint, Acrylic , "SourcePointLogical");
		SourcePointLogical->SetVisAttributes(PETcolour); //raphael(changed SourcePoint material from PVC to Acrylic)
		new G4PVPlacement(0, G4ThreeVector(0,0,0) , SourcePointLogical,"SourcePhysical",SourceFrameLogical, false, 0);

		// by bjoern
		// constructing source holder top plate
		if (ActivateHADESSourceHolder == true) {
		  
		  G4double SHLidInnerRadius = pow(pow(11.08,2)+pow(23.48,2),0.5)/2.*mm + 0.1*mm; // diagonal of source frame rectangular + epsilon 
                                                                                                 // raphael: changed slightly the values to the actual source dimensions
		  G4Tubs* SHLid = new G4Tubs("SHLid", SHLidInnerRadius, SHInnerRadius,
					     SHLidHeight / 2., 0, 360. * deg);
		  G4LogicalVolume *SHLidLogical = new G4LogicalVolume(SHLid, Acrylic,
								      "SHLidLogical");
		  SHLidLogical->SetVisAttributes(Acryliccolour);
		  new G4PVPlacement(0, *SourceVector + G4ThreeVector(0,0,fSourceThickness/2.-SHLidHeight/2.), SHLidLogical,
				    "SHLidPhysical", experimentalHallLogical, false, 0);
		}
		//end by bjoern
	}

	//#################################################
	//LPRIBox: Information by Erica (same as CERCA?????)
	else if(fSourceType=="Cerca")
	{
		MGLog(trace)<<"Loading "<<fSourceType<<" source... "<<endlog;

		/*
		 source thickness according to specs:
		 mass per area: 28 mg*cm-2
		 density PE (wiki): 0.915-0.97 g/cm3  //??? check if its PE!!!
		 thickness: 0.29-0.31 mm
		 average thickness: 0.3mm

		 The source is constructed with an outer Plastic ring that holds two layers of thin plastic (PE???) foil that are sandwitching a "very thin compact-grain 3mm layer of source material"

		 */
		G4double fSourceThickness = 0.3 *mm; // thickness of one side of foil
		G4double fSourceRingThickness = 3 *mm;// thickness of one side of foil

		//source location
		G4ThreeVector* SourceVector = new G4ThreeVector(0,0,CryostatHeight/2. + fSourceRingThickness/2. + SourceDistance);

		G4RotationMatrix* SourceRotation = new G4RotationMatrix();

		//source ring
		G4Tubs* SourceRing = new G4Tubs("SourceRing",
				16.001/2. * mm,
				25./2. * mm,
				fSourceRingThickness/2.,
				0.0 * deg,
				360 * deg);
		G4LogicalVolume* SourceRingLogical
		= new G4LogicalVolume(SourceRing, PE, "SourceRingLogical");
		new G4PVPlacement(SourceRotation, *SourceVector, SourceRingLogical,"SourceRingPhysical",experimentalHallLogical, false, 0);

		// source foil
		G4Tubs* SourceFoil = new G4Tubs("SourceFoil",
				0. * mm,
				16./2. * mm,
				fSourceThickness,
				0.0 * deg,
				360 * deg);
		G4LogicalVolume* SourceFoilLogical
		= new G4LogicalVolume(SourceFoil, PE, "SourceFoilLogical");
		new G4PVPlacement(0, *SourceVector , SourceFoilLogical,"SourceFoilPhysical",experimentalHallLogical, false, 0);

		// source point
		G4Tubs* SourcePoint = new G4Tubs("SourcePoint",
				0. * mm,
				3./2. * mm,
				0.005 * mm,
				0.0 * deg,
				360 * deg);
		G4LogicalVolume* SourcePointLogical
		= new G4LogicalVolume(SourcePoint, PE , "SourcePointLogical");
		SourcePointLogical->SetVisAttributes(PETcolour);
		new G4PVPlacement(0, G4ThreeVector(0,0,0) , SourcePointLogical,"SourcePhysical",SourceFoilLogical, false, 0);

	}

	//#################################################
	//HS7, HS8 source: Information from HEROICA wiki
	else if((fSourceType=="HS7")||(fSourceType=="HS7like"))
	{
		MGLog(trace)<<"Loading "<<fSourceType<<" source... "<<endlog;

		/*
		 source thickness according to specs:
		 mass per area: 28 mg*cm-2
		 density PE (wiki): 0.915-0.97 g/cm3  //??? check if its PE!!!
		 thickness: 0.29-0.31 mm
		 average thickness: 0.3mm

		 The source is constructed with an outer Plastic ring that holds two layers of thin plastic (PE???) foil that are sandwitching a "very thin compact-grain 3mm layer of source material"

		 */
		G4double fSourceThickness = 3. *mm; // thickness of one side of foil

		//source location
		G4ThreeVector* SourceVector = new G4ThreeVector(0,0,CryostatHeight/2. + fSourceThickness/2. + SourceDistance);

		G4RotationMatrix* SourceRotation = new G4RotationMatrix();

		//source coating
		G4Tubs* SourceCoating = new G4Tubs("SourceCoating",
				0. * mm,
				25./2. * mm,
				fSourceThickness/2.,
				0.0 * deg,
				360 * deg);

		G4LogicalVolume* SourceCoatingLogical
		= new G4LogicalVolume(SourceCoating, ssteel , "SourceCoatingLogical");
		new G4PVPlacement(SourceRotation, *SourceVector, SourceCoatingLogical,"SourceCoatingPhysical",experimentalHallLogical, false, 0);

		// source point
		G4Box *SourcePoint = new G4Box("SourcePoint",
				1./2. * mm,
				1./2. * mm,
				1./2. * mm);

		G4NistManager* manager = G4NistManager::Instance();
		G4Material* ceramic = manager->FindOrBuildMaterial("G4_SILICON_DIOXIDE");


		G4LogicalVolume* SourcePointLogical
		= new G4LogicalVolume(SourcePoint, ceramic, "SourcePointLogical");
		SourcePointLogical->SetVisAttributes(Ceramiccolour);
		new G4PVPlacement(0, G4ThreeVector(0,0,0) , SourcePointLogical,"SourcePhysical",SourceCoatingLogical, false, 0);


	}



	//#################################################
	//LNGSWireSource: Information by Nuno / Stefan
	else if((fSourceType=="LNGSWireSource")||(fSourceType=="HS6"))
	{
		MGLog(trace)<<"Loading "<<fSourceType<<" source... "<<endlog;
		// FIXME: This is not being used.
		// TODO:
		G4double SourceCoatingLength = 3.5*mm;
		G4double SourceCoatingRadius = 2.0*mm;

		G4ThreeVector* SourceVector = new G4ThreeVector(0,0,CryostatHeight/2. + SourceCoatingRadius + SourceDistance);

		G4RotationMatrix* SourceRotation = new G4RotationMatrix();
		SourceRotation ->rotateY(90*deg);

		// The source is composed of different regions (see schematic)
		// An inner cylinder of Epoxy
		// Another inner cylinder of ceramic
		// Covered by a shell of stainless steel

		// Build the source itself. Forget about the encasing.

		MGLog(warning) << "=========== Starting material loop :" << endlog;
		for (unsigned int im = 0; im < G4Material::GetNumberOfMaterials(); ++im ) {
			MGLog(warning) << "Material : " << (*(G4Material::GetMaterialTable()))[im]->GetName() << endlog;
		}
		MGLog(warning) << "=========== Finished material loop." << endlog;

		G4Tubs* CalibSourceTubs = new G4Tubs("SourceCaseTubs",
				0.0,
				1.0 * mm,
				3.5 * mm,
				0.0 * deg,
				360 * deg);

		G4LogicalVolume* steelSourceLogical =
		new G4LogicalVolume(CalibSourceTubs,G4Material::GetMaterial("Steel"),"SourceCaseLogisteelSourceLogicalcal");

		G4VPhysicalVolume* steelShieldingPhysical =
		new G4PVPlacement(SourceRotation, *SourceVector, steelSourceLogical,"SourceCase",experimentalHallLogical, false, 0);

		// Now build the internal parts of the source:
		//      The 1mm ceramic Th source
		//       The epoxy layer
		G4Tubs* CeramicSourceTubs = new G4Tubs("CalibSourceTubs",
				0.0,
				0.5 * mm,
				0.5 * mm,
				0.0 * deg,
				360 * deg);

		G4NistManager* manager = G4NistManager::Instance();
		G4Material* mat = manager->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

		G4LogicalVolume* ceramicSourceLogical =
		new G4LogicalVolume(CeramicSourceTubs,mat,"CalibSourceLogical");
		ceramicSourceLogical->SetVisAttributes(PETcolour);

		G4VPhysicalVolume* ceramicShieldingPhysical =
		new G4PVPlacement(0,G4ThreeVector(0,0,0) ,ceramicSourceLogical,"SourcePhysical",steelSourceLogical,false,0);

		// Add the epoxy layer to the botom of the source
		G4Tubs* epoxySourceTubs = new G4Tubs("SourceEpoxyTubs",
				0.0,
				0.5 * mm,
				1.5 * mm,
				0.0 * deg,
				360 * deg);

		G4LogicalVolume* epoxySourceLogical =
		new G4LogicalVolume(epoxySourceTubs,G4Material::GetMaterial("Epoxy_mod"),"SourceEpoxyLogical");

		G4VPhysicalVolume* epoxyShieldingPhysical =
		new G4PVPlacement(0,G4ThreeVector(0,0,-0.5*mm - 1.5*mm) ,epoxySourceLogical,"EpoxySource",steelSourceLogical,false,0);

	}

	//#################################################
	//Systematic error study for HS21 Am241 source geometry:
	else if((fSourceType=="SSG0")||(fSourceType=="SSG1")||(fSourceType=="SSG2")||(fSourceType=="SSG3")||(fSourceType=="SSG4")||(fSourceType=="SSG5")||(fSourceType=="SSG6"))
	{
		/*
		 SSG0: identical to HS21 (thickness=2.3*mm, r=0.5*mm)
		 SSG1: increase source volume (1.0*mm)
		 SSG2: decrease source volume (0.1*mm)
		 SSG3: increase source thickness (2.8*mm)
		 SSG4: decrease source thickness (1.8*mm)
		 SSG5: SourcePoint 1.0mm thickness 2.8*mm
		 SSG6: SourcePoint 0.1mm thickness 1.8*mm

		 */

		MGLog(trace)<<"Loading "<<fSourceType<<" source... "<<endlog;

		G4double fSourceThickness = 2.3*mm; //
		G4double fSourcePointRadius = 1./2. * mm;//

		if(fSourceType=="SSG1") {fSourcePointRadius = 1.0*mm;}
		else if(fSourceType=="SSG2") {fSourcePointRadius = 0.1*mm;}
		else if(fSourceType=="SSG3") {fSourceThickness = 2.8*mm;}
		else if(fSourceType=="SSG4") {fSourceThickness = 1.8*mm;}
		else if(fSourceType=="SSG5") {fSourcePointRadius = 1.0*mm; fSourceThickness = 2.8*mm;}
		else if(fSourceType=="SSG6") {fSourcePointRadius = 0.1*mm; fSourceThickness = 1.8*mm;}

		MGLog(trace)<<"fSourceThickness set to "<<fSourceThickness/1*mm<<"mm"<<endlog;
		MGLog(trace)<<"fSourcePointRadius set to "<<fSourcePointRadius/1*mm<<"mm"<<endlog;

		//source location
		G4ThreeVector* SourceVector = new G4ThreeVector(0,0,CryostatHeight/2. + fSourceThickness/2. + SourceDistance);

		G4RotationMatrix* SourceRotation = new G4RotationMatrix();

		//source frame
		G4Box *SourceFrame = new G4Box("SourceFrame",
				10./2. * mm,
				23./2. * mm,
				fSourceThickness/2.);

		G4LogicalVolume* SourceFrameLogical
		= new G4LogicalVolume(SourceFrame, PVC, "SourceFrameLogical");
		new G4PVPlacement(SourceRotation, *SourceVector, SourceFrameLogical,"SourceFramePhysical",experimentalHallLogical, false, 0);

		// source point
		G4Sphere* SourcePoint = new G4Sphere("SourcePoint",
				0. * mm,
				fSourcePointRadius,
				0.0 * deg,
				360 * deg,
				0.0 * deg,
				180 * deg);

		G4LogicalVolume* SourcePointLogical
		= new G4LogicalVolume(SourcePoint, PVC , "SourcePointLogical");
		SourcePointLogical->SetVisAttributes(PETcolour);
		new G4PVPlacement(0, G4ThreeVector(0,0,0) , SourcePointLogical,"SourcePhysical",SourceFrameLogical, false, 0);

	}


	else if(fSourceType=="flat")
	{
		MGLog(trace)<<"Loading "<<fSourceType<<" source... "<<endlog;

		/*
		- flat source with 0.1mm thickness of air with the diameter exactly matching the bege diameter
		- used for generating homogeneous photon filed
		 */
		G4double fSourceThickness = 0.1 *mm; // thickness of one side of foil

		//source location
		G4ThreeVector* SourceVector = new G4ThreeVector(0,0,CryostatHeight/2. + fSourceThickness/2. + SourceDistance);


		// source point
		G4Tubs* SourcePoint = new G4Tubs("SourcePoint",
				0. * mm,
				XtalDiameter/2.,
				fSourceThickness,
				0.0 * deg,
				360 * deg);
		G4LogicalVolume* SourcePointLogical
		= new G4LogicalVolume(SourcePoint, Air , "SourcePointLogical");
		SourcePointLogical->SetVisAttributes(PETcolour);
		new G4PVPlacement(0, *SourceVector , SourcePointLogical,"SourcePhysical",experimentalHallLogical, false, 0);

	}


	// ####################################################################
	// ------------------- lead castle ------------------------------------

	G4double cryostatDepth = 9.2 * cm;

	// no lead castle
	if (HADESLeadCastleType == 0) {

	}

	// lead castle without cu
	if ((HADESLeadCastleType == 1) || (HADESLeadCastleType == 2)) {
		// top of cryostat is at 1/3 of castle height
		G4ThreeVector* CastleVector = new G4ThreeVector(0, 0,
				0 + CryostatHeight / 2. + 30. / 2. * cm - cryostatDepth);

		//outer box
		G4Box *leadBox = new G4Box("leadBox", 30. / 2. * cm, 30. / 2. * cm,
				30. / 2. * cm);
		//inner cutout
		G4Box *leadBoxNeg = new G4Box("leadBoxNeg", (30. - 10.) / 2. * cm,
				(30. - 10.) / 2. * cm, (30. / 2. + 0.1) * cm);

		G4VSolid* leadBox_2 = new G4SubtractionSolid("leadBox_2", leadBox,
				leadBoxNeg, 0, G4ThreeVector(0, 0, 0));

		G4LogicalVolume* leadBoxLogical = new G4LogicalVolume(leadBox_2, Pb,
				"leadBoxLogical");
		leadBoxLogical->SetVisAttributes(Pbcolour);

		new G4PVPlacement(0, *CastleVector, leadBoxLogical, "leadBoxPhysical",
				experimentalHallLogical, false, 0);

		//top plate 3.5mm steel
		if (HADESLeadCastleType == 1) {

			G4Box *topLid = new G4Box("topLid", 30. / 2. * cm, 30. / 2. * cm,
					3.5 / 2. * mm);

			G4LogicalVolume* topLidLogical = new G4LogicalVolume(topLid, ssteel,
					"topLidLogical");
			topLidLogical->SetVisAttributes(Ssteelcolour);

			new G4PVPlacement(0,
					G4ThreeVector(0, 0,
							0 + CryostatHeight / 2. + 30. * cm - cryostatDepth
									+ 3.5 / 2. * mm), topLidLogical,
					"topLidPhysical", experimentalHallLogical, false, 0);

		}

	}

	// lead castle with cu inlet
	if (HADESLeadCastleType == 2) {
		// top of cryostat is at 1/3 of castle height
		G4ThreeVector* CastleVector = new G4ThreeVector(0, 0,
				0 + CryostatHeight / 2. + 30. / 2. * cm - 9.2 * cm);

		//outer box
		G4Box *CuBox = new G4Box("leadBox", (20. / 2.) * cm, (20. / 2.) * cm,
				(30. / 2.) * cm);
		//inner cutout
		G4Box *CuBoxNeg = new G4Box("leadBoxNeg", (20. - 6.) / 2. * cm,
				(20. - 6.) / 2. * cm, (30. / 2. + 0.1) * cm);

		G4VSolid* CuBox_2 = new G4SubtractionSolid("CuBox_2", CuBox, CuBoxNeg,
				0, G4ThreeVector(0, 0, 0));

		G4LogicalVolume* CuBoxLogical = new G4LogicalVolume(CuBox_2, Cu,
				"CuBoxLogical");
		CuBoxLogical->SetVisAttributes(Cucolour);

		new G4PVPlacement(0, *CastleVector, CuBoxLogical, "CuBoxPhysical",
				experimentalHallLogical, false, 0);

		//top plate 10mm Cu
		G4Box *topLid = new G4Box("topLid", 30. / 2. * cm, 30. / 2. * cm,
				10. / 2. * mm);

		G4LogicalVolume* topLidLogical = new G4LogicalVolume(topLid, Cu,
				"topLidLogical");
		topLidLogical->SetVisAttributes(Cucolour);

		new G4PVPlacement(0,
				G4ThreeVector(0, 0,
						0 + CryostatHeight / 2. + 30. * cm - cryostatDepth
								+ 10 / 2. * mm), topLidLogical,
				"topLidPhysical", experimentalHallLogical, false, 0);

	}

	// ##################################################################################################
	// ------------------- BRADY collimated Am241 source environment ------------------------------------
	if (fUseBRADYEnv) {

		if ((fUseCollimator == 1) || (fUseCollimator == 2)) {

			// BARDY Am Collimator Properties
			G4double collimatorLength = 65. * mm;
			G4double collimatorWidth = 30. * mm;
			G4double collimatorThickness = 30. * mm;

			G4double collimatorBeamRadius = 0.5 * mm;
			G4double collimatorBeamLength = 25.6 * mm;

			G4double collimatorY = 0., collimatorZ = 0.;

			// Cylindrical encapsulation volume
			G4double encapsRadius = 1. * mm;
			G4double encapsLength = 10. * mm;

			// Activity
			G4double sourceRadius = 0.5 * mm;
			G4double sourceLength = 2. * mm;
			G4double sourceZ = sourceLength / 2. - encapsLength / 2. + 0.2 * mm;

			// Collimator Volume
			G4Box* collimator = new G4Box("collimator", collimatorWidth / 2.,
					collimatorThickness / 2., collimatorLength / 2.);
			G4LogicalVolume* logicalCollimator = new G4LogicalVolume(collimator,
					Cu, "logicCollimator");
			logicalCollimator->SetVisAttributes(Cucolour);

			// Collimator Beam Volume
			G4Tubs* collimatorBeam = new G4Tubs("collimatorBeam", 0,
					collimatorBeamRadius, collimatorBeamLength / 2., 0, twopi);
			G4LogicalVolume* logicalCollimatorBeam = new G4LogicalVolume(
					collimatorBeam, Air, "logicCollimatorBeam");
			logicalCollimatorBeam->SetVisAttributes(Vacuumcolor);

			// Encapsulation Volume
			G4Tubs* encapsulation = new G4Tubs("encapsulation", 0, encapsRadius,
					encapsLength / 2., 0, twopi);
			G4LogicalVolume* logicalEncapsulation = new G4LogicalVolume(
					encapsulation, ssteel, "logicEncapsulation");
			logicalEncapsulation->SetVisAttributes(Ssteelcolour);

			// Source Volume
			G4Tubs* sourceVolume = new G4Tubs("sourceVolume", 0, sourceRadius,
					sourceLength / 2., 0, twopi);
			G4LogicalVolume* logicalSourceVolume = new G4LogicalVolume(
					sourceVolume, AmericiumOxide, "logicSourceVolume");
			logicalSourceVolume->SetVisAttributes(AmOxcolour);

			// Rotation Matrix
			G4RotationMatrix* rotateCollimator = new G4RotationMatrix();

			if (fUseCollimator == 1) //Side Scan
					{
				MGLog(trace)<< "Constructing Collimator for Side Scan: " << fCollimatorPosition/mm << " mm " << endlog;

				// Rotation Matrix
				rotateCollimator->rotateX(90.*deg);

				collimatorY = CryostatDiameter/2. + fCollimatorDistance + collimatorLength/2.;
				collimatorZ = CryostatHeight/2. - fCollimatorPosition;
			}
			else //Top Scan
			{
				MGLog(trace) << "Constructing Collimator for Top Scan: " << fCollimatorPosition/mm << " mm " << endlog;

				// Rotation Matrix
				rotateCollimator->rotateX(0.);

				collimatorY = fCollimatorPosition;
				collimatorZ = CryostatHeight/2. + fCollimatorDistance + collimatorLength/2.;
			}

			// Collimator Placement
			new G4PVPlacement(rotateCollimator,
					G4ThreeVector(0, collimatorY, collimatorZ),
					logicalCollimator, "physicalCollimator",
					experimentalHallLogical, false, 0);
			// Collimator Beam Placement
			new G4PVPlacement(0,
					G4ThreeVector(0, 0,
							collimatorBeamLength / 2. - collimatorLength / 2.),
					logicalCollimatorBeam, "physicalCollimatorBeam",
					logicalCollimator, false, 0);
			// Encapsulation Volume Placement
			new G4PVPlacement(0,
					G4ThreeVector(0, 0,
							encapsLength / 2. + collimatorBeamLength
									- collimatorLength / 2.),
					logicalEncapsulation, "physicalEncapsulation",
					logicalCollimator, false, 0);
			// Source Volume Placement
			new G4PVPlacement(0, G4ThreeVector(0, 0, sourceZ),
					logicalSourceVolume, "BRADYSourcePhysical",
					logicalEncapsulation, false, 0);

		} else if (fUseCollimator == 0) {
			MGLog(trace)<< "Collimator set to false." << endlog;
		}
		else {
			MGLog(trace) << "Unknown option for Collimator: " << fUseCollimator << endlog;
			MGLog(trace) << "Use: 0=false, 1=SideScan, 2=TopScan" << endlog;
		}

	}
	else {
		MGLog(trace) << "BRADY Environment set to false" << endlog;
	}

}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

G4float GEHADESBEGeTests::GetDistToUpperSurf(G4ThreeVector vec) {
	//  cryostat window ends with cryostat

	G4double zPositionBEGeUpperSurf = CryostatHeight / 2
			- CryostatWindowThickness - XtalDistanceToWindow;
	return zPositionBEGeUpperSurf - vec.getZ();
}

G4float GEHADESBEGeTests::GetDistToLowerSurf(G4ThreeVector vec) {
	G4double zPositionBEGeLowerSurf = CryostatHeight / 2
			- CryostatWindowThickness - XtalDistanceToWindow - XtalHeight;
	G4double hitRadius = pow(pow(vec.getX(), 2) + pow(vec.getY(), 2), 0.5);
	G4double retVal = 0 * cm;

	// decide if hit is within inner electrode outside. return physical distance only if outside
	if (hitRadius < XtalDitchOuterRadius) {
		retVal = -1 * cm;
	} else {
		retVal = vec.getZ() - zPositionBEGeLowerSurf;
	}

	return retVal;
}

G4float GEHADESBEGeTests::GetDistToSideSurf(G4ThreeVector vec) {
	G4double hitRadius = pow(pow(vec.getX(), 2) + pow(vec.getY(), 2), 0.5);
	return XtalDiameter / 2. - hitRadius;
}

G4float GEHADESBEGeTests::GetDistToCornerSurf(G4ThreeVector vec) {
	//strategy:
	// fist: calculate radius at heigth of hit position (C)
	// second: calculate distance to side surface at height of hit position (K)
	// third: calculate angle of corner (phi)
	// forth: calculate shortest distance to corner surface (DCS)

	G4double DCS; //Distance to Corner Surface

	if ((XtalCornerDiameter == 0) && (XtalCornerHeight == 0)) {
		DCS = -1 * cm;
	} else {
		//first:
		G4double C = XtalCornerDiameter / 2.
				+ (XtalDiameter / 2. - XtalCornerDiameter / 2.)
						/ XtalCornerHeight * GetDistToUpperSurf(vec);

		//second:
		G4double K = GetDistToSideSurf(vec) - (XtalDiameter / 2. - C);

		//third:
		G4double tanPhi = XtalCornerHeight
				/ (XtalDiameter / 2. - XtalCornerDiameter / 2.);

		//forth:
		DCS = sin(atan(tanPhi)) * K;
	}

	return DCS;
}

G4float GEHADESBEGeTests::GetDistToDitchSurf(G4ThreeVector vec) {
	//Improve me
	return 0;
}
