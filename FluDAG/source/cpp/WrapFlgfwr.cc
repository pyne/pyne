///////////////////////////////////////////////////////////////////
//
// WrapFlgfwr.cc - Paola Sala
//
// Wrapper for setting of fluka geometry flag
//
// created  12/Nov/09
//
//////////////////////////////////////////////////////////////////

#include "DagWrappers.hh"
// #include "FGeometryInit.hh"

void flgfwr ( int& flkflg )
{
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "=======FLGFWR =======" << std::endl;
#endif
/*
   static FGeometryInit * ptrGeoInit = FGeometryInit::GetInstance();
#ifdef DAGGEOMETRY_DEBUG
  int geoflg =   ptrGeoInit-> getGeoFlag();
  std::cout << "=======FLGFWR =======" << geoflg << std::endl;    
#endif
  ptrGeoInit->setGeoFlag(flkflg);//  geoflg = flkflg;
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "=======FLGFWR =======" << std::endl;
  geoflg =   ptrGeoInit-> getGeoFlag();
  std::cout << "=======FLGFWR =======" << geoflg << std::endl;    
#endif
*/
}
