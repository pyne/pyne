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

double posx;
double posy;
double posz;
double xyz[] = {posx, posy, posz};

double dirx;
double diry;
double dirz;

double lastRetStep;
int lastKnownRegion = -1;
int curReg, newReg, nextReg, originReg;;

char *getfileptr(std::string name);

int look( double& posx, double& posy, double& posz, double* dir, int& region);
void normal (double& posx, double& posy, double& posz, double *norml, int& regionForSense);

TEST(nrml_Test, spheres_sense)
{
   double retStep;
   posx = 0.0;
   posy = 0.0;
   posz = 0.0;
   
   dirx = 1.0;
   diry = 0.0;
   dirz = 0.0;
   double dir[] = {dirx, diry, dirz};
   double* norml = new double(3); 
   int flagErr;
   
   // Nested bricks: faces at 5, 7.5, and 10
   // Concentric spheres at 12, 15, 20
   char *fileptr = getfileptr("../cases/spheres.h5m");
   std::cout << "********* Testing Concentric Spheres ************" << std::endl;

   // Load the file and create the dag structure
   cpp_dagmcinit(fileptr,0,1);
   
   xyz = {posx, posy, posz}; 
   // Find out what volume we are in: lastKnownRegion doesn't matter unless we are on a boundary
   // curRegion should be 1 for nested cubes, 2 for nested spheres
   curReg = look(posx, posy, posz, dir, lastKnownRegion );
   originReg = curReg;

   // Fire a ray from the given point in the given direction
   // calls ray_fire; sets retStep, nextReg and global next_surf
   g1_fire(curReg, xyz, dir, retStep, nextReg); 
   // nrmlwr(posx, posy, posz, dirx, diry, dirz, norml, curReg, newReg, flagErr);
   normal(posx, posy, posz, norml, curReg);
   // Potential Tests: 
   //   o verify retStep = {12,5} at origin for spheres and bricks
   //   o verify curReg = {2,1} at origin for spheres and bricks. 
   //   o Reset the position to just a little off the origin, within the normaliztion amount
   //     (1e-16) and verify nothing changes

   // Set file global for later use
   // lastRetStep = retStep;
   if (true)
   {
      std::cout<<"============= Result of g1_fire =============="<<std::endl;
      std::cout << "Position " << posx << " " << posy << " " << posz << std::endl;
      std::cout << "Point is in region = " << curReg<< std::endl;
      std::cout << "Direction vector " << dirx << " " << diry << " " << dirz << std::endl;
      std::cout << "retStep = " << retStep << std::endl;
      std::cout << "Point heading to region = " << nextReg << std::endl;
   }
   std::cout << std::endl;


   // Testing calls (not otherwise necessary) 
   // Test:  confirm near origin we are still in the volume of the origin 
   std::cout << "Testing position close to origin " << std::endl;
   posx += 1.0e-2;
   lastKnownRegion = curReg;  // Copy position of previous region
   curReg = look(posx, posy, posz, dir, lastKnownRegion);
   EXPECT_EQ(originReg, curReg);

   // Position is near the origin.  Nothing should change for the bricks.  For the sphere,
   // it appears that slight changes in position pick up different facets on which to 
   // calculate the normal.
   normal(posx, posy, posz, norml, curReg);
   // Get the sense of the normal with respect to the current region
   int sense = getSense(curReg);
   std::cout << "Sense of next_surf with respect to the point is " << sense << std::endl;
   std::cout << std::endl;
   EXPECT_EQ(1, sense);

   std::cout << "Fire the ray in the opposite direction.  " << std::endl;
   dir[0] *= -1.0;
   xyz = {posx, posy, posz}; 
   g1_fire(curReg, xyz, dir, retStep, nextReg);
   if (true)
   {
      std::cout<<"============= Result of g1_fire =============="<<std::endl;
      std::cout << "Position " << posx << " " << posy << " " << posz << std::endl;
      std::cout << "Point is in region = " << curReg<< std::endl;
      std::cout << "Direction vector " << dir[0] << " " << dir[1] << " " << dir[2] << std::endl;
      std::cout << "retStep = " << retStep << std::endl;
      std::cout << "Point heading to region = " << nextReg << std::endl;
   }
   std::cout << std::endl;
   normal(posx, posy, posz, norml, curReg);
   sense = getSense(curReg);
   std::cout << "Sense of next_surf with respect to the point is " << sense << std::endl;
}

TEST(nrml_test, test2)
{ 
   double* norml = new double(3); 
   double retStep;
   int flagErr;
   double dir[] = {-dirx, diry, dirz};
   // Add distance to surface to current (already offset) point should put us into the next region
   std::cout << std::endl;
   posx = lastRetStep + .1;
   lastKnownRegion = curReg;  // Copy position of previous region
   curReg = look(posx, posy, posz, dir, lastKnownRegion);
   // Expect curReg to = previous nextReg from knowledge of last g1_fire, but maybe not because we are such
   // a short distance beyond the boundary.
   // calls ray_fire, sets retStep, nextReg and global next_surf
   // double original_retStep = retStep;  // save retstep
   g1_fire(curReg, xyz, dir, retStep, nextReg); 
   if (true)
   {
      std::cout<<"============= Result of g1_fire =============="<<std::endl;
      std::cout << "Position " << posx << " " << posy << " " << posz << std::endl;
      std::cout << "Direction vector " << dirx << " " << diry << " " << dirz << std::endl;
      std::cout << "Point is in region = " << curReg<< std::endl;
      std::cout << "retStep = " << retStep << std::endl;
      std::cout << "Point heading to region = " << nextReg << std::endl;
   }
   nrmlwr(posx, posy, posz, dirx, diry, dirz, norml, curReg, newReg, flagErr);
   int sense = getSense(curReg);
   std::cout << "Sense of next_surf with respect to the point is " << sense << std::endl;
   std::cout << std::endl;

   /////////////////////////////////////////////////////////////////////////////   
   // Move this to another test
   /*
   std::cout << std::endl;
   posx += retStep;
   lastKnownRegion = curReg;  // Copy position of previous region
   curReg = look(posx, posy, posz, dir, lastKnownRegion);
   // Expect curReg to = previous nextReg from knowledge of last g1_fire, but maybe not because we are such
   // a short distance beyond the boundary.
   // calls ray_fire, sets retStep, nextReg and global next_surf
   g1_fire(curReg, xyz, dir, retStep, nextReg); 
   if (true)
   {
      std::cout<<"============= Result of g1_fire =============="<<std::endl;
      std::cout << "Position " << posx << " " << posy << " " << posz << std::endl;
      std::cout << "Direction vector " << dirx << " " << diry << " " << dirz << std::endl;
      std::cout << "Point is in region = " << curReg<< std::endl;
      std::cout << "retStep = " << retStep << std::endl;
      std::cout << "Point heading to region = " << nextReg << std::endl;
   }
   nrmlwr(posx, posy, posz, dirx, diry, dirz, norml, curReg, newReg, flagErr);
   */
}

//---------------------------------------------------------------------------//
// getfileptr
//---------------------------------------------------------------------------//
// helper function
char *getfileptr(std::string name)
{
   char* fileptr = new char(name.length()+1);
   std:strcpy(fileptr, name.c_str());
}

//---------------------------------------------------------------------------//
// look
//---------------------------------------------------------------------------//
// Local wrapper for fortran-called, fluka-set lkwr
// This function signature shows what parameters are being used in our wrapper implementation
// ASSUMES:  position is not on a boundary
int look( double& posx, double& posy, double& posz, double* dir, int& region)
{
   int flagErr;
   int lattice_dummy;  // not used
   int nextRegion;     // looked at iff on boundry, but not set.
   lkwr(posx, posy, posz, dir, region, lattice_dummy, nextRegion, flagErr, lattice_dummy);
   return nextRegion;
}

//---------------------------------------------------------------------------//
// normal
//---------------------------------------------------------------------------//
// Local wrapper for fortran-called, fluka-set nrmlwr
// This function signature shows what parameters are being used in our wrapper implementation
// ASSUMES:  no ray history
// Notes
// - direction is not taken into account 
// - curRegion is not currently used.  It is expected to be implemented as a check
//   on what the sign of the normal should be.  It is used in a call to DAG->surface_sense
void normal (double& posx, double& posy, double& posz, double *norml, int& curRegion)
{
   int dummyFlagErr, dummyReg;
   double dummyDirx, dummyDiry, dummyDirz;
   nrmlwr(posx, posy, posz, dummyDirx, dummyDiry, dummyDirz, norml, curRegion, dummyReg, dummyFlagErr);
}
