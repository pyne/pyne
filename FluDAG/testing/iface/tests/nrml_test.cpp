//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   /testing/iface/test/nrml_test.cpp
 * \author Julie C. Zachman
 * \date   Wed Apr 3 
 * \brief  Unit Tests for Issue
 * \note   Test
 */
//---------------------------------------------------------------------------//
// $Id:
//---------------------------------------------------------------------------//
#include "gtest/gtest.h"
#include "DagWrappers.hh"
#include "fluka_funcs.h"
#include "MBInterface.hpp"
#include "DagMC.hpp"
#include "moab/Types.hpp"

using moab::DagMC;

double pSx;
double pSy;
double pSz;

double pVx;
double pVy;
double pVz;
int oldReg;

char *getfileptr(std::string name);

TEST(nrml_Test, RayFire)
{
   pSx = 0;
   pSy = 0;
   pSz = 0;
   double xyz[] = {pSx, pSy, pSz};
   
   pVx = 1;
   pVy = 0;
   pVz = 0;
   double dir[] = {pVx, pVy, pVz};
   
   int oldReg = 2;
   double retStep;
   int newReg = 0;

   char *fileptr;
   // fileptr = getfileptr("../iface/cases/test.h5m");
   fileptr = getfileptr("../iface/cases/test.h5m");

   std::string fname = "../iface/cases/test.h5m";

   // Load the file and create the dag structure
   // cpp_dagmcinit(fileptr,0,1);
   cpp_dagmcinit(fileptr,0,1);
   // Fire a ray from the given point in the given direction
   // calls ray_fire, sets retStep and global next_surf
   
   g1_fire(oldReg, xyz, dir, retStep, newReg); 
   if (true)
   {
      std::cout<<"============= Result of g1_fire =============="<<std::endl;
      std::cout << "Position " << pSx << " " << pSy << " " << pSz << std::endl;
      std::cout << "Direction vector " << pVx << " " << pVy << " " << pVz << std::endl;
      std::cout << "Oldreg = " << oldReg << std::endl;
      std::cout << "retStep = " << retStep << std::endl;
      std::cout << "Newreg = " << newReg << std::endl;
   }
   // Now get the normal of next_surf, as set by the g1_fire call
   // Question:  do I need to pretend fluka updated the particle to the boundary?
   //            In this case I should reset pSx, pSy, pSz to the next_surf
   double* norml = new double(3); 
   int flagErr;
   nrmlwr(pSx, pSy, pSz, pVx, pVy, pVz, norml, 0,0, flagErr);

   std::cout<<"============= Norml after g1_fire =============="<<std::endl;
   std::cout << "Normal: " << norml[0] << ", " << norml[1] << ", " << norml[2]  << std::endl;
   std::cout << "retStep = " << retStep << std::endl;

   pSx += retStep + .001;
   nrmlwr(pSx, pSy, pSz, pVx, pVy, pVz, norml, 0,0, flagErr);

   std::cout<<"============= Norml after g1_fire =============="<<std::endl;
   std::cout << "Normal: " << norml[0] << ", " << norml[1] << ", " << norml[2]  << std::endl;
   std::cout << "retStep = " << retStep << std::endl;

   EXPECT_TRUE(true);
}

char *getfileptr(std::string name)
{
   char* fileptr = new char(name.length()+1);
   std:strcpy(fileptr, name.c_str());
   return fileptr;
}
/*
//---------------------------------------------------------------------------//
// nrmlwr(..)
//---------------------------------------------------------------------------//
/// From Flugg Wrappers WrapNorml.cc
void nrmlwr(double& pSx, double& pSy, double& pSz,
            double& pVx, double& pVy, double& pVz,
	    double* norml, const int& oldReg, 
	    const int& newReg, int& flagErr)
{
  if(debug)
    {
      std::cout << "============ NRMLWR-DBG =============" << std::endl;
    }

  //dummy variables
  flagErr=0;
  double abc[] = {pSx, pSy, pSz};
  double xyz[3]; //tmp storage of position
  //MBEntityHandle surf = 0;
  xyz[0]=pSx,xyz[1]=pSy,xyz[2]=pSz;
  MBErrorCode ErrorCode = DAG->get_angle(next_surf,xyz,norml); // get the angle
  if(ErrorCode != MB_SUCCESS)
    {
      std::cout << "Could not determine normal" << std::endl;
      flagErr = 2;
      return;
    }

  //return normal:
  //norml[0]=pVx;
  //norml[1]=pVy;
  //norml[2]=pVz;

  // to test:  create a simple model:  dag-> ray_fire in order to get next_sruf.  
// then call normlwr with next_surf and it should return an opposite-pointing vector
// ON the surface, normlwr should components of rnorml should be 0 or near
// PAST the surface, norml components should point AWAY from current position
  
  if(debug)
    {
      std::cout << "Normal: " << norml[0] << ", " << norml[1] << ", " << norml[2]  << std::endl;
      std::cout << "out of nrmlwr " << std::endl;
    }
  
  return;
}
*/
