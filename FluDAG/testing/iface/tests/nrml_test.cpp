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

// #define DAG DagMC::instance()

double posx_orig;
double posx;
double posy;
double posz;

double dirx;
double diry;
double dirz;
int oldReg;

char *getfileptr(std::string name);
/*
MBEntityHandle next_srf;
void nrml(double& pSx, double& pSy, double& pSz,
            double& dirx, double& diry, double& dirz,
            double* norml, const int& oldReg,
            const int& newReg, int& flagErr);

void g1_rayfire(int& oldRegion, double point[], double dir[], double& retStep,  int& newRegion);
*/
void delta_to_surface(double posx, double posy, double posz,
                      double dirx, double diry, double dirz,
                      double delta, int& oldReg, double retStep, int newReg,
                      int& flagErr, bool print_all);

TEST(nrml_Test, RayFire)
{
   // double delta = 0.00;
   posx_orig = 0.0;
   posx = posx_orig;
   posy = 0.0;
   posz = 0.0;
   double xyz[] = {posx, posy, posz};
   
   dirx = 1.0;
   diry = 0.0;
   dirz = 0.0;
   double dir[] = {dirx, diry, dirz};
   
   int oldReg = 2;
   double retStep;
   int newReg;
   int flagErr;
   int lattice_dummy;

   char *fileptr;
   fileptr = getfileptr("../cases/test.h5m");

   // Load the file and create the dag structure
   cpp_dagmcinit(fileptr,0,1);
   lkwr(posx, posy, posz, dir, oldReg, lattice_dummy, newReg, flagErr, lattice_dummy);
   
   // Fire a ray from the given point in the given direction
   // calls ray_fire, sets retStep and global next_surf
   // std::cout << next_surf << std::endl;
   g1_fire(oldReg, xyz, dir, retStep, newReg); 
   // std::cout << next_surf << std::endl;

   if (true)
   {
      std::cout<<"============= Result of g1_fire =============="<<std::endl;
      std::cout << "Position " << posx << " " << posy << " " << posz << std::endl;
      std::cout << "Direction vector " << dirx << " " << diry << " " << dirz << std::endl;
      std::cout << "Oldreg = " << oldReg << std::endl;
      std::cout << "retStep = " << retStep << std::endl;
      std::cout << "Newreg = " << newReg << std::endl;
   }
   std::cout << std::endl;

   double delta = 1.0e-16;
   delta_to_surface(posx, posy, posz, dirx, diry, dirz, delta, \
                    oldReg, retStep, newReg, flagErr, false);

   delta = 1.0e-15;
   delta_to_surface(posx, posy, posz, dirx, diry, dirz, delta, \
                    oldReg, retStep, newReg, flagErr, false);

   delta = 1.0e-14;
   delta_to_surface(posx, posy, posz, dirx, diry, dirz, delta, \
                    oldReg, retStep, newReg, flagErr, false);

   delta = 1.0e-13;
   delta_to_surface(posx, posy, posz, dirx, diry, dirz, delta, \
                    oldReg, retStep, newReg, flagErr, false);

   delta = 1.0e-12;
   delta_to_surface(posx, posy, posz, dirx, diry, dirz, delta, \
                    oldReg, retStep, newReg, flagErr, false);

   delta = 1.0e-11;
   delta_to_surface(posx, posy, posz, dirx, diry, dirz, delta, \
                    oldReg, retStep, newReg, flagErr, false);

   delta = 1.0e-10;
   delta_to_surface(posx, posy, posz, dirx, diry, dirz, delta, \
                    oldReg, retStep, newReg, flagErr, false);

   delta = 1.0e-9;
   delta_to_surface(posx, posy, posz, dirx, diry, dirz, delta, \
                    oldReg, retStep, newReg, flagErr, false);

   delta = 1.0e-8;
   delta_to_surface(posx, posy, posz, dirx, diry, dirz, delta, \
                    oldReg, retStep, newReg, flagErr, false);

   delta = 1.0e-7;
   delta_to_surface(posx, posy, posz, dirx, diry, dirz, delta, \
                    oldReg, retStep, newReg, flagErr, false);

   delta = 1.0e-6;
   delta_to_surface(posx, posy, posz, dirx, diry, dirz, delta, \
                    oldReg, retStep, newReg, flagErr, false);

   delta = 1.0e-5;
   delta_to_surface(posx, posy, posz, dirx, diry, dirz, delta, \
                    oldReg, retStep, newReg, flagErr, false);

   delta = 1.0e-4;
   delta_to_surface(posx, posy, posz, dirx, diry, dirz, delta, \
                    oldReg, retStep, newReg, flagErr, false);

   delta = 1.0e-3;
   delta_to_surface(posx, posy, posz, dirx, diry, dirz, delta, \
                    oldReg, retStep, newReg, flagErr, false);

   EXPECT_TRUE(true);
}

char *getfileptr(std::string name)
{
   char* fileptr = new char(name.length()+1);
   std:strcpy(fileptr, name.c_str());
   return fileptr;
}


void delta_to_surface(double posx, double posy, double posz,
                      double dirx, double diry, double dirz,
                      double delta, int& oldReg, double retStep, int newReg,
                      int& flagErr, bool print_all)
{
   // Now get the normal of next_surf, as set by the g1_fire call
   // Question:  do I need to pretend fluka updated the particle to the boundary?
   //            In this case I should reset posx, posy, posz to the next_surf
   double* norml = new double(3); 
   std::cout << "For delta = " << delta << ": Position still at origin " << posx << " " << posy << " " << posz << std::endl;
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

/*
void nrml(double& pSx, double& pSy, double& pSz,
            double& dirx, double& diry, double& dirz,
            double* norml, const int& oldReg,
            const int& newReg, int& flagErr)
{
  if(true)
    {
      std::cout << "============ NRMLWR-DBG =============" << std::endl;
    }

  //dummy variables
  flagErr=0;
  double xyz[3]; //tmp storage of position
  //MBEntityHandle surf = 0;
  xyz[0]=pSx,xyz[1]=pSy,xyz[2]=pSz;
  MBErrorCode ErrorCode = DAG->get_angle(next_srf,xyz,norml); // get the angle
  if(ErrorCode != MB_SUCCESS)
    {
      std::cout << "Could not determine normal" << std::endl;
      flagErr = 2;
      return;
    }

  //return normal:
  //norml[0]=dirx;
  //norml[1]=diry;
  //norml[2]=dirz;

  // to test:  create a simple model:  dag-> ray_fire in order to get next_sruf.  
// then call normlwr with next_surf and it should return an opposite-pointing vector
// ON the surface, normlwr should components of rnorml should be 0 or near
// PAST the surface, norml components should point AWAY from current position

  if(true)
    {
      std::cout << "Normal: " << norml[0] << ", " << norml[1] << ", " << norml[2]  << std::endl;
      std::cout << "out of nrmlwr " << std::endl;
    }

  return;
}
///////                 End nrmlwr
/////////////////////////////////////////////////////////////////////

//---------------------------------------------------------------------------//
// g1(int& old Region, int& newRegion)
//---------------------------------------------------------------------------//
void g1_rayfire(int& oldRegion, double point[], double dir[], double& retStep,  int& newRegion)
{
  MBEntityHandle vol = DAG->entity_by_index(3,oldRegion);

  double next_surf_dist;
  MBEntityHandle newvol = 0;

  MBErrorCode result = DAG->ray_fire(vol, point, dir, next_srf, next_surf_dist );
  retStep = next_surf_dist;

  MBErrorCode rval = DAG->next_vol(next_srf,vol,newvol);

  newRegion = DAG->index_by_handle(newvol);
  return;
}
///////                 End g1wr and g1
/////////////////////////////////////////////////////////////////////
*/

/*
//---------------------------------------------------------------------------//
// nrmlwr(..)
//---------------------------------------------------------------------------//
/// From Flugg Wrappers WrapNorml.cc
void nrmlwr(double& posx, double& posy, double& posz,
            double& dirx, double& diry, double& dirz,
	    double* norml, const int& oldReg, 
	    const int& newReg, int& flagErr)
{
  if(debug)
    {
      std::cout << "============ NRMLWR-DBG =============" << std::endl;
    }

  //dummy variables
  flagErr=0;
  double abc[] = {posx, posy, posz};
  double xyz[3]; //tmp storage of position
  //MBEntityHandle surf = 0;
  xyz[0]=posx,xyz[1]=posy,xyz[2]=posz;
  MBErrorCode ErrorCode = DAG->get_angle(next_surf,xyz,norml); // get the angle
  if(ErrorCode != MB_SUCCESS)
    {
      std::cout << "Could not determine normal" << std::endl;
      flagErr = 2;
      return;
    }

  //return normal:
  //norml[0]=dirx;
  //norml[1]=diry;
  //norml[2]=dirz;

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
