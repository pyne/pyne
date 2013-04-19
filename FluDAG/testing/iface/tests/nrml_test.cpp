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

TEST(nrml_Test, sense)
{
   double retStep;
   posx = 0.0;
   posy = 0.0;
   posz = 0.0;
   
   dirx = 1.0;
   diry = 0.0;
   dirz = 0.0;
   double dir[] = {dirx, diry, dirz};
   
   // Concentric bricks: faces at 5, 7.5, and 10
   char *fileptr = getfileptr("../cases/test.h5m");

   // Load the file and create the dag structure
   cpp_dagmcinit(fileptr,0,1);
   
   // Find out what volume we are in: lastKnownRegion doesn't matter unless we are on a boundary
   // curRegion should be 1 for nested cubes, 2 for nested spheres
   curReg = look(posx, posy, posz, dir, lastKnownRegion );
   originReg = curReg;

   // Fire a ray from the given point in the given direction
   // calls ray_fire, sets retStep, newReg and global next_surf
   g1_fire(curReg, xyz, dir, retStep, nextReg); 
   // set global for later use
   lastRetStep = retStep;
   if (true)
   {
      std::cout<<"============= Result of g1_fire =============="<<std::endl;
      std::cout << "Position " << posx << " " << posy << " " << posz << std::endl;
      std::cout << "Direction vector " << dirx << " " << diry << " " << dirz << std::endl;
      std::cout << "Point is in region = " << curReg<< std::endl;
      std::cout << "retStep = " << retStep << std::endl;
      std::cout << "Point heading to region = " << nextReg << std::endl;
   }
   std::cout << std::endl;

   // Reset the position to just a little off the origin, within the normaliztion amount
   double* norml = new double(3); 
   posx += 1.0e-2;
   std::cout << "Position close to origin " << posx << " " << posy << " " << posz << std::endl;

   // Testing calls (not otherwise necessary: confirm we are still in the volume of the origin 
   lastKnownRegion = curReg;  // Copy position of previous region
   curReg = look(posx, posy, posz, dir, lastKnownRegion);
   EXPECT_EQ(originReg, curReg);

   // Position is near the origin.  Nothing should change
   int flagErr;
   nrmlwr(posx, posy, posz, dirx, diry, dirz, norml, curReg, newReg, flagErr);
   // Get the sense of the normal with respect to the current region
   int sense = getSense(curReg);
      std::cout << "Sense of next_surf with respect to the point is " << sense << std::endl;
   EXPECT_EQ(1, sense);

   dir[0] *= -1.0;
   g1_fire(curReg, xyz, dir, retStep, nextReg);
   nrmlwr(posx, posy, posz, dirx, diry, dirz, norml, curReg, newReg, flagErr);
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
int look( double& posx, double& posy, double& posz, double* dir, int& region)
{
   int flagErr;
   int lattice_dummy;  // not used
   int nextRegion;     // looked if on boundry, but not set.
   lkwr(posx, posy, posz, dir, region, lattice_dummy, nextRegion, flagErr, lattice_dummy);
   return nextRegion;
}
//---------------------------------------------------------------------------//
// delta_to_surface
//---------------------------------------------------------------------------//
/*
void delta_to_surface(double posx, double posy, double posz,
                      double dirx, double diry, double dirz,
                      double delta, int& oldReg, double retStep, int newReg,
                      int& flagErr, bool print_all)
{
   // Now get the normal of next_surf, as set by the g1_fire call
   // Question:  do I need to pretend fluka updated the particle to the boundary?
   //            In this case I should reset posx, posy, posz to the next_surf
   double* norml = new double(3); 
   std::cout << "Position still at origin " << posx << " " << posy << " " << posz << std::endl;
   // Position is still at the origin
   nrmlwr(posx, posy, posz, dirx, diry, dirz, norml, oldReg, newReg, flagErr);
   
   if (print_all)
   {
      std::cout << "Oldreg = " << oldReg << std::endl;
      std::cout << "Newreg = " << newReg << std::endl;
      std::cout << std::endl;
   }
   // Save, but not used
   double nrml_orig_pos[] =  {norml[0], norml[1], norml[2]};

   // Update position to delta before surfac3
   std::cout<<"============= update position to " << delta << \
              "  before surface =============="<<std::endl;
   posx = posx_orig + (retStep - delta);
   nrmlwr(posx, posy, posz, dirx, diry, dirz, norml, oldReg, newReg, flagErr);
   if (print_all)
   {
   std::cout << "Position " << posx << " " << posy << " " << posz << std::endl;
      std::cout << "Oldreg = " << oldReg << std::endl;
      std::cout << "Newreg = " << newReg << std::endl;
      std::cout << std::endl;
   }
   double nrml_before_pos[] =  {norml[0], norml[1], norml[2]};  // save for comparing

   // Update position to the surface
   std::cout<<"============= update position at surface =============="<<std::endl;
   posx = posx_orig + retStep;
   nrmlwr(posx, posy, posz, dirx, diry, dirz, norml, oldReg, newReg, flagErr);
   if (print_all)
   {
   std::cout << "Position " << posx << " " << posy << " " << posz << std::endl;
      std::cout << "Oldreg = " << oldReg << std::endl;
      std::cout << "retStep = " << retStep << std::endl;
      std::cout << "Newreg = " << newReg << std::endl;
      std::cout << std::endl;
   }
   double nrml_at_pos[] =  {norml[0], norml[1], norml[2]};    // save for comparing

   // Update position to delta after surface
   std::cout<<"============= update position to " << delta << \
              " beyond surface =============="<<std::endl;
   posx = posx_orig + (retStep + delta);
   nrmlwr(posx, posy, posz, dirx, diry, dirz, norml, oldReg, newReg, flagErr);
   if (print_all)
   {
   std::cout << "Position " << posx << " " << posy << " " << posz << std::endl;
      std::cout << "Oldreg = " << oldReg << std::endl;
      std::cout << "retStep = " << retStep << std::endl;
      std::cout << "Newreg = " << newReg << std::endl;
      std::cout << std::endl;
   }
   double nrml_after_pos[] =  {norml[0], norml[1], norml[2]};  // save for comparing

   std::string sense1;
   if (nrml_at_pos[0] * nrml_before_pos[0] >= 0)
   {
       sense1 = "same";
   }
   else
   {
       sense1 = "OPPOSITE";
   }

   std::cout << "Normal AT surface -  Normal BEFORE surface: " << std::endl;
   std::cout << (nrml_at_pos[0] - nrml_before_pos[0]) << ", " << \
                (nrml_at_pos[1] - nrml_before_pos[1]) << ", " << \
                (nrml_at_pos[2] - nrml_before_pos[2]) << " in " << \
                sense1 << " direction." << std::endl;

   std::string sense2;
   if (nrml_at_pos[0] * nrml_after_pos[0] >= 0)
   {
       sense2 = "same";
   }
   else
   {
       sense2 = "OPPOSITE";
   }
   std::cout << "Normal AT surface -  Normal AFTER surface: " << std::endl;
   std::cout << (nrml_at_pos[0] - nrml_after_pos[0]) << ", " << \
                (nrml_at_pos[1] - nrml_after_pos[1]) << ", " << \
                (nrml_at_pos[2] - nrml_after_pos[2]) << " in " << \
                sense2 << " direction." << std::endl;
   std::cout << std::endl;
}
*/
