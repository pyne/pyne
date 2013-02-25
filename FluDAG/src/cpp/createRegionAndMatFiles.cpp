#include <iostream>
#include <stdlib.h>
//#include <string>
#include <stdio.h>
#include <fstream>

#include "DagWrapUtils.hh"

using namespace moab;

#define DAG DagMC::instance()


//*****************************************************************************


void createFlukaMatFile() 
{
  
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "==> FluDAG createFlukaMatFile()" << std::endl;
  std::cout << "================== FILEWR =================" << std::endl;
#endif 


  std::ofstream vos("Volumes_index.inp");
  PrintEntityRegionNames(vos);
  vos.close();

  //Materials and compounds
  BuildMaterialTables();
  std::ofstream fos("flukaMat.inp");  
/*
  PrintMaterialTables(fos);
  PrintAssignmat(fos);
  PrintMagneticField(fos);
  fos.close();
*/

#ifdef DAGGEOMETRY_DEBUG
  std::cout << "<== FluDAG createFlukaMatFile()" << std::endl;
#endif
}

////////////////////////////////////////////////////////////////////////
//  			makeRegionName()
//////////////////////////////////////////////////////////////////////
///////				Create a name for Entity l
////// No Maps!
/////////////////////////////////////////////////////////////////////
std::string makeRegionName(int l)
{

  std::string VVname;
  std::string Vname;
  EntityHandle entity = NULL;
  entity = DAG->entity_by_id(3, l);
  // Don't add 1
  int iRegion = l;
  std::cout << iRegion << " index,name " << std::endl;
  char vname[8];
  sprintf(vname,"%-8u",iRegion);
  Vname.replace(0,8,vname);
  std::cout<<iRegion<<" vname" << vname <<" Vname " << Vname<<std::endl;
  unsigned int found=Vname.rfind(" ",7);
  // take out blanks
  while (found<Vname.size()-1)
  {
      found=Vname.rfind(" ",7);
      if (found<Vname.size()-1)
      {
	  std::string temp=Vname.substr(found+1);
	  Vname.replace(found, temp.size(), temp);
	  Vname.erase(Vname.size()-1);
	  std::cout << Vname << std::endl;
      } 
  }

    unsigned int nameSize=Vname.size();
    if (nameSize > 8 ) {VVname=Vname.substr(0,7);}
    else {VVname=Vname;}
    std::cout << "VVname "<<VVname<< std::endl;
    // ToDo:  Re-implement the maps if we need to check that the name is unique and there
    // isn't another way.
    // jcz - ask about this
    // check that name is unique, if not add numbering
/*
    unsigned int matrep = 1;
    unsigned int ii=VVname.size()+3;
    unsigned int  newSize = ii < 8 ? ii : 8;
    bool old=true;
    char smatrep[3];
       while (old)
       {
	 old=false;
         for ( DRegionIterator i=dRegionNameMap.begin(); i!=dRegionNameMap.end();i++)
         {
           sprintf(smatrep,"%03u",matrep);
           Vname=(*i).second;
           if (Vname==VVname)
           {
   	      old=true;
   	      if (VVname.size()<=5)
              {
   	          VVname+=smatrep;
              }
   	      else
              {
   	          VVname.resize(newSize);
   	          VVname.replace(VVname.size()-3, 3,smatrep);
              }
	   matrep++;
           }
         }
       }
      std::cout<< "VVname "<< VVname<< std::endl;
*/
      return VVname;
}
//////////////////////////////////////////////////////////////////////
///////			End makeRegionName
/////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
///////			PrintEntityRegionNames()
/////////////////////////////////////////////////////////////////////
void PrintEntityRegionNames(std::ostream& os)
{
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "==> FluDAG PrintEntityNames()" << std::endl;
#endif
  PrintHeader(os, "VOLUMES, and Fake Names");
  std::string Vname;
  unsigned int numVol = DAG->num_entities(3);
  for(unsigned int l=0; l < numVol; l++) 
  {
    int iRegion = l+1;
    Vname = makeRegionName(iRegion);
    //Print index and region name in some fixed format
    writeRegionLine(os, iRegion, Vname);
  }
  int iRegion = numVol + 1;
  Vname = "BLACKHOL";
  writeRegionLine(os, iRegion, Vname);

#ifdef DAGGEOMETRY_DEBUG
  std::cout << "<== FluDAG PrintEntityRegionNames()" << std::endl;
#endif
}

/////////////////////////////////////////////////////////////////////
// Write one line of the *.inp file
// Convenience function that is re-used
/////////////////////////////////////////////////////////////////////
void writeRegionLine(std::ostream& os, int iRegion, std::string name)
{
    os.setf(std::ios::left, std::ios::adjustfield);
    os << setw10 << iRegion;
    // jcz - Money call:  how to get the name out of MBEntity?
    // ToDo:  replace "entname" if possible
    // os << std::setw(20) << entity->GetName() << std::setw(20) << "";
    os << std::setw(20) << "entname"  << std::setw(20) << "";
    os << std::setw(5) << name << std::setw(5) << "";
    // ToDo:  duplicate this for DAG calls, if possible and necessary
    //If volume is a replica... print some more stuff
/*
    if(ptrVol->IsReplicated()) 
    {
      EAxis axis;
      int nRep = -1;
      G4double width = -1;
      G4double offset = -1;
      G4bool consum = false;
      ptrVol->GetReplicationData(axis, nRep, width, offset, consum);
      os.setf(std::ios::left, std::ios::adjustfield);
      os << setw10 << "Repetion Nb: " << std::setw(3) << nRep;
    }
*/
    os << std::endl;
}
//////////////////////////////////////////////////////////////////////
///////			End PrintEntityRegionNames()
/////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
// 
// ToDo:  implement this if necessary.
/*
int FGeometryInit::GetRegionFromName(const char* volName) const {
  for (RegionIterator i = fRegionVolumeMap.begin(); 
       i != fRegionVolumeMap.end(); 
       i++) {
    
    //Get info in the map
    G4VPhysicalVolume* ptrVol = (*i).first;
    if (ptrVol->GetName() == volName)
      return ((*i).second);
  }
  return -1;
}
*/



////////////////////////////////////////////////////////////////////////
// 
void BuildMaterialTables() 
{
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "==> FluDAG BuildMaterialTables()" << std::endl;
#endif

  //some terminal printout also
  std::cout << "\t* Storing information..." << std::endl;

  //The logic is the folloing:
  //Get the Material Table and:
  // 1) For materials with density <= 1.00e-10*g/cm3 assign vacuum
  // 2) For each single element material build a material equivalent
  // 3) For the rest:
  //   3.a) Build materials for each not already known element
  //   3.b) Build the compound out of them

  // initialize with predefined materials
  initmat = InitFlukaMat();
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "end init predef mat  "  << std::endl;
#endif
  //Get the Material Table and iterate
  const G4MaterialTable* matTable = DagG4Material::GetMaterialTable();
  for (MatTableIterator i = matTable->begin(); i != matTable->end(); i++) {

    //Get some basic material information
    G4Material* material = (*i);
    DString matName = material->GetName();
    const G4double matDensity = material->GetDensity();
    const int nMatElements  = material->GetNumberOfElements();
#ifdef DAGGEOMETRY_DEBUG
    std::cout << " treating material " << matName    << std::endl;
#endif

    std::cout << " mat " << matName 
	   << ": dens. = " << matDensity/(g/cm3) << "g/cm3"
	   << ", nElem = " << nMatElements << std::endl;

    // 1) For materials with density <= 1.00e-10*g/cm3 assign vacuum
    //    FlukaMaterial* is  in that case
    if (matDensity <= 1.00e-10*g/cm3) {
#ifdef DAGGEOMETRY_DEBUG
    std::cout << " vacuum?  "<< matDensity    << std::endl;
#endif
      DString elemName("VACUUM");
      FlukaMaterial *flukamat = FlukaMaterial::GetFlukaMaterial(elemName);
      G4FlukaMaterialMap[material] = flukamat;
      std::cout << "\t\t  Stored as " << flukamat->GetRealName() << std::endl;
    }
    // 2) For each single element material build a material equivalent
    else if (nMatElements == 1) {
#ifdef DAGGEOMETRY_DEBUG
    std::cout << " single element "     << std::endl;
#endif
      
      FlukaMaterial *flukamat = 
	BuildFlukaMaterialFromElement(material->GetElement(0),
				      matDensity);
      
      G4FlukaMaterialMap[material] = flukamat;
      std::cout << "  Stored as " << flukamat->GetRealName() << std::endl;
      
    } //else if (material->GetNumberOfElements() == 1)
    
    // 3) For the rest:
    //   3.a) Build materials for each not already known element
    //   3.b) Build the compound out of them
    else {
#ifdef DAGGEOMETRY_DEBUG
      std::cout << " not  vacuum : call Comp. "<< matDensity/(g/cm3)   << std::endl;
#endif
      FlukaCompound* flukacomp = 
	BuildFlukaCompoundFromMaterial(material);
      G4FlukaCompoundMap[material] = flukacomp;
      std::cout << "\t\t  Stored as " << flukacomp->GetRealName() << std::endl;
    } //else for case 3)
  } //for (materials)
  
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "<== Flugg FGeometryInit::BuildMaterialTables()" << std::endl;
#endif
}


/*
FlukaMaterial* 
FGeometryInit::BuildFlukaMaterialFromElement(const G4Element* element,
					     G4double matDensity) {
#ifdef DAGGEOMETRY_DEBUG
  cout << "==> Flugg FGeometryInit::BuildFlukaMaterialFromElement()" 
	 << endl;
#endif

  //Get element and its properties
  DString elemName(ToFlukaString(element->GetName()));
      FlukaMaterial* flukamat = FlukaMaterial::GetFlukaMaterial(elemName);
  //Check for isotopes
      int nIsotopes = element->GetNumberOfIsotopes();
#ifdef DAGGEOMETRY_DEBUG
      cout << " elem : " << elemName << " nIsot " << nIsotopes 	
	     << " density"  << matDensity << endl;

#endif
    if (nIsotopes == 0) {

  #ifdef DAGGEOMETRY_DEBUG
      cout << " NIso=0"<< endl;
#endif
      // matDensity =0  means element in compound. 
//  if (matDensity != 0 || (matDensity == 0 && flukamat == 0)) {
      G4bool lnewmat = true;
      if ( flukamat != 0 ) 
	{
       if (matDensity == 0 ) {
	 lnewmat =  false ; }
       else {
         G4double de =  flukamat->GetDensity();
	 // if previous dens=0, and new one differs, keep material 
         if ( de != 0 ){ 
#ifdef DAGGEOMETRY_DEBUG
	 DString fluna = flukamat->GetRealName();
	  G4double dG = matDensity/(g/cm3);
	  cout << " alredy there?  " << endl;
	cout << " fluka name " << fluna << " fluka d. " <<de << " G4 de  " << dG << endl;
#endif
         G4double dd = fabs (de - matDensity/(g/cm3) ) / de;
         if ( dd < 0.001 ) lnewmat = false;
         }
       } }
      // new material iso=o
       if (lnewmat)  
    {     
#ifdef DAGGEOMETRY_DEBUG
	 cout << " lnewmat = true " << endl;
#endif
      G4double elemA = element->GetA()/g;
      G4double elemZ = element->GetZ();
      
      if (elemA != int(elemA) && elemZ != int(elemZ)) {
	cout << "WARNING: Element \'" << elemName 
	       << "\' has non integer Z (" << elemZ << ") or A (" 
	       << elemA << ")"
	       << endl;
              }
    
        flukamat = new FlukaMaterial(elemName,
				   int(elemZ),
				   elemA,
				   matDensity/(g/cm3));
    }} //end nIso=0
    else if (nIsotopes == 1) {
  #ifdef DAGGEOMETRY_DEBUG
      cout << " NIso=1"<< endl;
   #endif
      const G4Isotope* isotope = element->GetIsotope(0);
      flukamat = BuildFlukaMaterialFromIsotope(isotope, matDensity);
    }
    else {
  #ifdef DAGGEOMETRY_DEBUG
      cout << " NIso>1"<< endl;
   #endif
      FlukaCompound *flucomp = BuildFlukaCompoundFromElement(element,
							     matDensity);
      flukamat = flucomp->GetFlukaMaterial();      
    }
  

  return flukamat;
  
#ifdef DAGGEOMETRY_DEBUG
  cout << "<== Flugg FGeometryInit::BuildFlukaMaterialFromElement()" 
	 << endl;
#endif
}
*/


/*
FlukaMaterial* 
FGeometryInit::BuildFlukaMaterialFromIsotope(const G4Isotope* isotope,
					     G4double matDensity) {
#ifdef DAGGEOMETRY_DEBUG
  cout << "==> Flugg FGeometryInit::BuildFlukaMaterialFromIsotope()" 
	 << endl;
#endif
  DString isoName(ToFlukaString(isotope->GetName()));
  FlukaMaterial* flukamat = FlukaMaterial::GetFlukaMaterial(isoName);
  int isoN = isotope->GetN();
  int isoF = 0;
  if (flukamat != 0)   isoF = flukamat->GetN(); 
  //  if (matDensity != 0 || (matDensity == 0 && flukamat == 0)) {
  if (flukamat == 0 || (isoF != isoN && flukamat != 0)) {
    G4double isoA = (isotope->GetA())/(g);
    int isoZ= isotope->GetZ();
    flukamat = new FlukaMaterial(isoName,
				 isoZ,
				 isoA,
				 matDensity/(g/cm3),
				 isoN);
  }

  return flukamat;

#ifdef DAGGEOMETRY_DEBUG
  cout << "==> Flugg FGeometryInit::BuildFlukaMaterialFromIsotope()" 
	 << endl;
#endif
}

FlukaCompound* 
FGeometryInit::BuildFlukaCompoundFromMaterial(const G4Material* material) {
#ifdef DAGGEOMETRY_DEBUG
  cout << "==> Flugg FGeometryInit::BuildFlukaCompoundFromMaterial()" 
	 << endl;
#endif
  //Material properties
  const G4double* elemFractions = material->GetFractionVector();
  const int nMatElements  = material->GetNumberOfElements();
  const G4double matDensity = material->GetDensity();
  DString matName(ToFlukaString(material->GetName()));
  FlukaCompound* flukacomp = new FlukaCompound(matName, matDensity/(g/cm3),
					       nMatElements);
  for (int i = 0; i < nMatElements; i++) {
#ifdef DAGGEOMETRY_DEBUG
    cout << "treating element n. " << i 
	 << endl;
    const G4Element* element = material->GetElement(i);
    DString elemName = element->GetName();
      DString uffa(ToFlukaString(elemName));
    cout << "its name " << elemName 
	   << " converted to " <<  uffa 
	 << endl;
#endif
    FlukaMaterial *flukamat = 
      BuildFlukaMaterialFromElement(material->GetElement(i), 0.0);
    
    flukacomp->AddElement(flukamat->GetIndex(), -elemFractions[i]);
    
  } //for (elements)

  return flukacomp;

#ifdef DAGGEOMETRY_DEBUG
  cout << "<== Flugg FGeometryInit::BuildFlukaCompoundFromMaterial()" 
	 << endl;
#endif
}
*/

/*
FlukaCompound* 
FGeometryInit::BuildFlukaCompoundFromElement(const G4Element* element,
					     G4double matDensity) {
#ifdef DAGGEOMETRY_DEBUG
  cout << "==> Flugg FGeometryInit::BuildFlukaCompoundFromElement()" 
	 << endl;
#endif
  int nIsotopes = element->GetNumberOfIsotopes();
  //fraction of nb of atomes per volume (= volume fraction?)
  const G4double* isoAbundance =  element->GetRelativeAbundanceVector();
  DString elemName(ToFlukaString(element->GetName()));

  //Material properties
  FlukaCompound* flukacomp = new FlukaCompound(elemName, matDensity/(g/cm3),
					       nIsotopes);
  for (int i = 0; i < nIsotopes; i++) {
    FlukaMaterial *flukamat = 
      BuildFlukaMaterialFromIsotope(element->GetIsotope(i), 0.0);
    
    flukacomp->AddElement(flukamat->GetIndex(), isoAbundance[i]);
    
  } //for (elements)

  return flukacomp;

#ifdef DAGGEOMETRY_DEBUG
  cout << "<== Flugg FGeometryInit::BuildFlukaCompoundFromElement()" 
	 << endl;
#endif
}
*/


int InitFlukaMat()
{
  int NumFlukaMat = 25 ;  
  DString FlukaNames[25] = { "BLCKHOLE" , "VACUUM", 
				 "HYDROGEN" ,   "HELIUM"  , 
                                  "BERYLLIU"  ,   "CARBON"  , 
	"NITROGEN"  ,   "OXYGEN"  ,   "MAGNESIU"  ,   "ALUMINUM"  , 
	  "IRON"  ,   "COPPER"  ,   "SILVER"  ,   "SILICON"  , 
	    "GOLD"  ,   "MERCURY"  ,   "LEAD"  ,   "TANTALUM"  ,
	      "SODIUM"  ,   "ARGON"  ,   "CALCIUM "  ,   "TIN"  ,
	"TUNGSTEN"  ,"TITANIUM"  ,"NICKEL" } ;
  double AFluka[25] = {        0.0     ,    0.0     ,
        1.00794 ,   4.002602,   9.012182,  12.0107  ,
       14.0067  ,  15.9994  ,  24.3050  ,  26.981538,
       55.845   ,  63.546   , 107.8682  ,  28.0855  ,
      196.96655 , 200.59    , 207.2     , 180.9479  ,
       22.989770,  39.948   ,  40.078   , 118.710   ,
       183.84    ,  47.867   ,  58.6934 };
  int ZFluka[25] = {       0, 0,    1 ,   2 ,   4 ,   6 ,
        7 ,   8 ,  12 ,  13 ,
       26 ,  29 ,  47 ,  14 ,
       79 ,  80 ,  82 ,  73 ,
       11 ,  18 ,  20 ,  50 ,
       74 ,  22 ,  28  };
  double RhoFluka[25] ={      0.0     ,    0.0     ,
       0.0837e-3,   0.166e-3,   1.848 ,   2.000  ,
        1.17e-3,   1.33e-3,   1.740 ,   2.699,
        7.874,   8.960,  10.500,   2.329,
       19.320,  13.546,  11.350,  16.654,
        0.971,   1.66e-3,   1.550,   7.310,
	19.300,   4.540,   8.902 };

  for ( int i=0; i != 25; i++ ) 
  {
    DString elemName = FlukaNames[i] ;
    double elemA = AFluka[i];
    int    elemZ = ZFluka[i];
    double matDensity = RhoFluka[i];
    FlukaMaterial *flukamat = new FlukaMaterial(elemName,
				   elemZ,
				   elemA,
				   matDensity);
     std::cout << " predef fluka mat " << flukamat->GetRealName() << std::endl;
  }

return  NumFlukaMat;
}


/*
void FGeometryInit::PrintMaterialTables(std::ostream& os) {
#ifdef DAGGEOMETRY_DEBUG
  cout << "==> Flugg FGeometryInit::PrintMaterialTables()" << endl;
#endif
  //Print Header
  PrintHeader(os, "GEANT4 MATERIALS AND COMPOUNDS");
  
  //And some more stuff  
  size_t nIsotopes = G4Isotope::GetNumberOfIsotopes();
  size_t nElements = G4Element::GetNumberOfElements();
  size_t nMaterials = G4Material::GetNumberOfMaterials();

  os << "* In Geant4 there are " << nMaterials << " materials" << endl;
  os << "* In Geant4 there are " << nElements  << " elements"  << endl;
  os << "* In Geant4 there are " << nIsotopes  << " isotopes"  << endl;

  //Materials
  cout << "\t* Printing FLUKA materials..." << endl;
  FlukaMaterial::PrintMaterialsByIndex(os);
  //FlukaMaterial::PrintMaterialsByName(os);

  //Compounds
  cout << "\t* Printing FLUKA compounds..." << endl;
  FlukaCompound::PrintCompounds(os);

#ifdef DAGGEOMETRY_DEBUG
  cout << "<== Flugg FGeometryInit::PrintMaterialTables()" << endl;
#endif
}
*/
////////////////////////////////////////////////////////////////////////
// 
/*
void PrintAssignmat(std::ostream& os) 
{
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "==> FluDAG PrintAssignmat()" << std::endl;
#endif

  //Find number of Volumes in physical volume store
  // G4PhysicalVolumeStore * pVolStore = G4PhysicalVolumeStore::GetInstance();
  // unsigned int numVol = pVolStore->size();
  
  unsigned int numVol = DAG->num_entities(3);

  std::cout << "\t* DAG has " << numVol << " volumes. " << std::endl;
  std::cout << "\t* Printing ASSIGNMAT..." << std::endl;

  EntityHandle entity;

  PrintHeader(os,"FluDAG MATERIAL ASSIGNMENTS");
  for(unsigned int l=0; l < numVol; l++) 
  {

    // Get each of the physical volumes
    //G4VPhysicalVolume * physicalvol = (*pVolStore)[l];
    entity = DAG->entity_by_id(3, l);

    // Get index for that volume
    //int iFlukaRegion = fRegionVolumeMap[physicalvol];
    int iRegion   = dRegionVolumeMap[&entity];
  

    // Get its name
    //DString Vname = fRegionNameMap[iFlukaRegion];
    std::string Vname = dRegionNameMap[iRegion];

    dRegionVolumeMap[&entity] = iRegion;

    //Find G4 material and navigate to its fluka compound/material
    // G4LogicalVolume * logicalVol = physicalvol->GetLogicalVolume();
    // G4Material* material = logicalVol->GetMaterial();
    // G4String G4volname = logicalVol->GetName();
    // jcz Don't have an answer to the above 4 lines:  maybe we don't need them

    //    cout << "l " << l << "volname " << G4volname << endl;
    int matIndex = 2;
    std::string mName ;

    if (G4FlukaCompoundMap[material]) 
    {
    matIndex = G4FlukaCompoundMap[material]->GetIndex();
     //   cout << "compound, index  " <<matIndex<< endl;
    mName =  G4FlukaCompoundMap[material]->GetRealName(); 
    }

    if (G4FlukaMaterialMap[material])
    {
    matIndex = G4FlukaMaterialMap[material]->GetIndex();
    // cout << "material, index  " <<matIndex<< endl;
    mName =  G4FlukaMaterialMap[material]->GetRealName(); 
    }
    //    cout << "mName " <<mName << endl;
//  jcz:  The following two if's assign matIndex, which is only used in a
          commented cout.
    if (G4FlukaCompoundMap[material]) 
    {
      matIndex = G4FlukaCompoundMap[material]->GetIndex();
    }
    if (G4FlukaMaterialMap[material])
    {
      matIndex = G4FlukaMaterialMap[material]->GetIndex();
    }

    //Find if there is a magnetic field in the region
    //check if Magnetic Field is present in the region

    G4double flagField = 0.0;
    G4FieldManager * pMagFieldMan = logicalVol->GetFieldManager();
    if(pMagFieldMan && pMagFieldMan->GetDetectorField())
      flagField = 1.0;

    //Print card
    os << setw10 << "ASSIGNMAT ";
    os.setf(static_cast<std::ios::fmtflags>(0),std::ios::floatfield);
    //os << setw10 << setfixed << G4double(matIndex);
       os << setw10  << mName;
//#ifdef NO_NAMES
    double iRegion_dble = (double)iRegion;
    os << setw10 << setfixed << iRegion_dble;
//#else
//    os << setw10 << setfixed << Vname;
//#endif
    os << setw10 << " ";
    // -- ahimmel add --
    os << setw10 << " ";
    // -- ahimmel end -- 
    // jcz - 'flagField' is associated with magnetic field
    double flagField = 0.0;
    os << setw10 << setfixed << flagField;
    os << std::endl;
  }
  //assign material 1 to black-hole=n.vol+1
  os << setw10 <<"ASSIGNMAT ";
  os << setw10 <<"1.0";
  os.setf(static_cast<std::ios::fmtflags>(0),std::ios::floatfield);
  os << setw10 << setfixed <<double(numVol+1);
  os << setw10 <<"0.0";
  os << setw10 <<"0.0";
  os << std::endl;



#ifdef DAGGEOMETRY_DEBUG
  std::cout << "==> FluDAG PrintAssignmat()" << std::endl;
#endif
}
*/
/*
void FGeometryInit::PrintMagneticField(std::ostream& os) {
#ifdef DAGGEOMETRY_DEBUG
  cout << "==> Flugg FGeometryInit::PrintMagneticField()" << endl;
#endif

  cout << "\t* Printing Magnetic Field..." << endl;

  if(fTransportationManager->GetFieldManager()->DoesFieldExist()) {
    
    //get magnetic field pointer
    const G4Field * pMagField = 
      fTransportationManager->GetFieldManager()->GetDetectorField();     
    
    
    if(pMagField) {
      G4double B[3];
      B[0] = 0.;
      B[1] = 0.;
      B[2] = 0.;
      //Check if it can be made a uniform magnetic field
      const G4UniformMagField *pUnifMagField = 
	dynamic_cast<const G4UniformMagField*>(pMagField);
      if(pUnifMagField) {
	G4double point[4]; //it is not really used
	pUnifMagField->GetFieldValue(point,B);
        B[0] = B[0]/tesla;
        B[1] = B[1]/tesla;
        B[2] = B[2]/tesla;
 
      }
      else {
	cout << "WARNING: No Uniform Magnetic Field found." << endl;
	cout << "         Manual intervention might be needed." << endl;
      }
      //write MGNFIELD card 
      PrintHeader(os,"GEANT4 MAGNETIC FIELD");
      os << setw10 << "MGNFIELD  ";
      os << setw10 << "";
      os << setw10 << "";
      os << setw10 << "";
      os.setf(static_cast<std::ios::fmtflags>(0),std::ios::floatfield);
      os << setw10 << setfixed
	 << std::setprecision(4) << B[0]
	 << setw10 << B[1]
	 << setw10 << B[2]
	 << endl;
    }
    else
      cout << "\t  No detector field found... " << endl;
  } // end if magnetic field
  else
    cout << "\t  No field found... " << endl;

#ifdef DAGGEOMETRY_DEBUG
  cout << "<== Flugg FGeometryInit::PrintMagneticField()" << endl;
#endif
}
*/

////////////////////////////////////////////////////////////////////////
// PrintHeader
////////////////////////////////////////////////////////////////////////
std::ostream& PrintHeader(std::ostream& os, const char* title) {
  os << "*\n" << "*\n" << "*\n";
  os << "*********************  " << title << " *********************\n"
     << "*\n";
  os << "*...+....1....+....2....+....3....+....4....+....5....+....6....+....7..."
     << std::endl;
  os << "*" << std::endl;

  return os;
}

