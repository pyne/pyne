
// FluDAG tag 

///////////////////////////////////////////////////////////////////
//
// WrapMag.hh - Sara Vanini
//
// Wrapper for geometry tracking in magnetic field: returns magnetic 
// field values in a given position.
//
// modified 26/X/1998
// modified 18/XI/1999
// modified 24.10.01: by I. Hrivnacova
//   functions declarations separated from implementation
//   (moved to Wrappers.hh);
//
/////////////////////////////////////////////////////////////////


#include "DagWrappers.hh"
// #include "FGeometryInit.hh"
// #include "globals.hh"
// -- ahimmel add --
// #include "G4VPhysicalVolume.hh"
// #include "G4PhysicalVolumeStore.hh"
// -- ahimmel end --

void fldwr(const double& pX, const double& pY, const double& pZ,
            double& cosBx, double& cosBy, double& cosBz, 
            double& Bmag, int& reg, int& idiscflag)

{
  //flag
#ifdef DAGGEOMETRY_DEBUG
  std::cout<<"================== MAGFLD ================="<<std::endl;
#endif 
/* 
  //Geoinit Pointer
  FGeometryInit * ptrGeoInit=FGeometryInit::GetInstance();
  
  //get FieldManager, Field pointers for magnetic field handling
  G4FieldManager * pFieldMgr = ptrGeoInit->getFieldManager();
  const G4Field * ptrField = pFieldMgr->GetDetectorField();
  
   // -- ahimmel add --
    if (!pFieldMgr->DoesFieldExist()) {
  G4PhysicalVolumeStore * pVolStore = G4PhysicalVolumeStore::GetInstance();
  G4VPhysicalVolume * ptrPhysReg = (*pVolStore)[reg-1];
  G4LogicalVolume* ptrLogReg = ptrPhysReg->GetLogicalVolume();
 
  // re-get FieldManager, Field pointers
  pFieldMgr = ptrLogReg->GetFieldManager();
  ptrField = pFieldMgr->GetDetectorField();
  } 
    // -- ahimmel end --
    
////compute field
  double point[3];
  point[0] = pX*10.;
  point[1] = pY*10.;
  point[2] = pZ*10.;
  double B[3];
  ptrField->GetFieldValue(point,B);
  double Bnor = sqrt(sqr(B[0])+sqr(B[1])+sqr(B[2]));
  if(Bnor) {
    cosBx = B[0]/Bnor;
    cosBy = B[1]/Bnor;
    cosBz = B[2]/Bnor;
  }
  else {
    cosBx = 0;
    cosBy = 0;
    cosBz = 1;
  }
  
  Bmag = Bnor/tesla;
  idiscflag = 0;
*/
}




