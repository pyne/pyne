//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   /testing/iface/test/nrml_spheres_test.cpp
 * \author Julie C. Zachman
 * \date   Wed Apr 3 
 * \brief  Unit Tests for concentric sphere geometry
 * \note   
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

double spheres_lastRetStep;

TEST(nrml_spheres_Test, sense)
{
   int lastKnownRegion = -1;
   int curReg, newReg, nextReg, originReg;;

   double propStep, retStep;
   // Position at zero
   double posx = 0.0;
   double posy = 0.0;
   double posz = 0.0;
   double xyz[] = {posx, posy, posz};
   
   // Direction along x
   double dirx = 1.0;
   double diry = 0.0;
   double dirz = 0.0;
   double dir[] = {dirx, diry, dirz};
   double* norml = new double(3); 
   int flagErr;
   
   // Concentric spheres at 12, 15, 20
   std::string name = "../cases/spheres.h5m";
   char* fileptr = new char(name.length()+1);
   std:strcpy(fileptr, name.c_str());
   std::cout << "********* Testing Concentric Spheres ************" << std::endl;

   // Load the file and create the dag structure
   cpp_dagmcinit(fileptr,0,1);
   
   // Find out what volume we are in: lastKnownRegion doesn't matter unless we are on a boundary
   // curRegion should be 1 for nested cubes, 2 for nested spheres
   curReg = look(posx, posy, posz, dir, lastKnownRegion );
   originReg = curReg;

   // Fire a ray from the given point in the given direction
   // calls ray_fire; sets retStep, nextReg and global next_surf
   g1_fire(curReg, xyz, dir, propStep, retStep, nextReg); 
   spheres_lastRetStep = retStep;  // save for later testing
   int error = normal(posx, posy, posz, norml, curReg);
   EXPECT_EQ(0,error);
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
      std::cout << "propStep = " << propStep << std::endl;
      std::cout << "retStep = " << retStep << std::endl;
      std::cout << "Point heading to region = " << nextReg << std::endl;
   }
   std::cout << std::endl;


   // Testing calls (not otherwise necessary) 
   // Test:  confirm near origin we are still in the volume of the origin 
   std::cout << "Testing position close to origin " << std::endl;
   posx += 1.0e-2;
   xyz[0] = posx;             // Update position vector
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
   g1_fire(curReg, xyz, dir, propStep, retStep, nextReg);
   if (true)
   {
      std::cout<<"============= Result of g1_fire =============="<<std::endl;
      std::cout << "Position " << posx << " " << posy << " " << posz << std::endl;
      std::cout << "Point is in region = " << curReg<< std::endl;
      std::cout << "Direction vector " << dir[0] << " " << dir[1] << " " << dir[2] << std::endl;
      std::cout << "propStep = " << propStep << std::endl;
      std::cout << "retStep = " << retStep << std::endl;
      std::cout << "Point heading to region = " << nextReg << std::endl;
   }
   std::cout << std::endl;
   normal(posx, posy, posz, norml, curReg);
   sense = getSense(curReg);
   std::cout << "Sense of next_surf with respect to the point is " << sense << std::endl;
}

TEST(nrml_spheres_Test, test2)
{ 
   int lastKnownRegion = -1;
   int curReg, newReg, nextReg, originReg;;
   double propStep, retStep;

   // Position at zero
   double posx = 0.0;
   double posy = 0.0;
   double posz = 0.0;
   double xyz[] = {posx, posy, posz};
   
   // Direction along -x
   double dirx = 1.0;
   double diry = 0.0;
   double dirz = 0.0;
   double dir[] = {-dirx, diry, dirz};
   double* norml = new double(3); 
   int flagErr;
   
   // Add distance to surface to current (already offset) point should put us into the next region
   std::cout << std::endl;
   posx = spheres_lastRetStep + .1;
   lastKnownRegion = curReg;  // Copy position of previous region
   curReg = look(posx, posy, posz, dir, lastKnownRegion);
   // Expect curReg to = previous nextReg from knowledge of last g1_fire, but maybe not because we are such
   // a short distance beyond the boundary.
   // calls ray_fire, sets retStep, nextReg and global next_surf
   // double original_retStep = retStep;  // save retstep
   g1_fire(curReg, xyz, dir, propStep, retStep, nextReg); 
   if (true)
   {
      std::cout<<"============= Result of g1_fire =============="<<std::endl;
      std::cout << "Position " << posx << " " << posy << " " << posz << std::endl;
      std::cout << "Direction vector " << dirx << " " << diry << " " << dirz << std::endl;
      std::cout << "Point is in region = " << curReg<< std::endl;
      std::cout << "propStep = " << propStep << std::endl;
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

