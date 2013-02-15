
///////////////////////////////////////////////////////////////////
//
// WrapInit.hh - Sara Vanini
//
// Wrapper for geometry initialisation.
//
// modified 12-IV-2000
// modified 24.10.01: by I. Hrivnacova
//   functions declarations separated from implementation
//   (moved to Wrappers.hh);
//
//////////////////////////////////////////////////////////////////

/*
#include <iostream>
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <fstream>
*/

#include "DagWrapUtils.hh"
#include "DagWrappers.hh"

using namespace moab;
#define DAG DagMC::instance()

void jomiwr(int & nge, const int& lin, const int& lou, int& flukaReg)
{
//flag
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "================== JOMIWR =================" << std::endl;
#endif 
  // return FLUGG code to fluka
  nge = 3;

  //Geoinit Pointer
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "\t *==> JOMIWR: Getting FGeometry..." << std::endl;
#endif 
//  FGeometryInit * ptrGeoInit=FGeometryInit::GetInstance();
  
  //initialize geometry:construct detector and set world volume
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "\t *==> JOMIWR: Setting the detector..." << std::endl;
#endif 
//  ptrGeoInit->setDetector();
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "\t *==> JOMIWR: Setting mother volume..." << std::endl;
#endif 
//  ptrGeoInit->setMotherVolume(); 
  
  //close geometry for optimization
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "\t *==> JOMIWR: Closing geometry..." << std::endl;
#endif 
//  ptrGeoInit->closeGeometry();
  
  //initialize wrappers utility histories at the beginning of run and set flag
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "\t *==> JOMIWR: InitHistories..." << std::endl;
#endif 
//  ptrGeoInit->InitHistories();
  
  //initialize lattice array
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "\t *==> JOMIWR: Init lattice array..." << std::endl;
#endif 
//  ptrGeoInit->InitJrLtGeantArray();
  
  //initialize debug-array
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "\t *==> JOMIWR: Init debug array..." << std::endl;
#endif 
//  ptrGeoInit->InitHistArray();
  
  //create Fluka material cards in flukaMat.inp file
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "\t *==> JOMIWR: Init fluka materials...STUBBED OUT FOR NOW jcz" << std::endl;
#endif 
//  ptrGeoInit->createFlukaMatFile();
//    createFlukaMatFile();
  
  //returns number of volumes + 1
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "\t *==> JOMIWR: Returning..." << std::endl;
#endif 
  // G4PhysicalVolumeStore * pVolStore = G4PhysicalVolumeStore::GetInstance();
  // int numVol = int(pVolStore->size());
  unsigned int numVol = DAG->num_entities(3);
  flukaReg = numVol + 1;
	
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "Number of volumes + 1: " << flukaReg << std::endl;
  std::cout << "================== Out of JOMIWR =================" << std::endl;
#endif
}



