// Define G4GEOMETRY_DEBUG for debugging information on cout

// #include "FGeometryInit.hh"
// #include "MyDetectorConstruction.hh"

#define flukam flukam_

/*
extern "C" int LOOKDB (double X, double Y, double Z, double MXERLM)
{
	return 1;
}
*/
extern "C" void flukam(const int &flag);

int main() {

//  FGeometryInit* theFGeometryInit = FGeometryInit::GetInstance();
  
// theFGeometryInit
//  ->setDetConstruction(new MyDetectorConstruction());

//flag for geometry:
// 1 for GEANT4
// 0 for FLUKA
// 2 for Rubia
// 3 for Dagmc
    const int flag = 3;

//call fortran
    flukam(flag);

//end
  return 0;
}




