//  DagSolid_test.cpp
#include <gtest/gtest.h>
#include <cmath>
#include <cassert>

#include <iostream>
#include "DagSolid.hh"
#include "G4TessellatedSolid.hh"

class DagSolidTest : public ::testing::Test
{
  protected:
};

/*
 * Point in volume test
 */
TEST_F(DagSolidTest,test_1) {
  DagMC* dagmc = DagMC::instance(); // create dag instance

  // dag_volumes
  const char* h5mfilename = "/data/opt/DAGMC/Geant4/test/geometry/test_geom.h5m";
  dagmc->load_file(h5mfilename,0);
  dagmc->init_OBBTree();

  // new volume
  DagSolid* vol_1 = new DagSolid("vol_1",dagmc,1);
  // sample position
  G4ThreeVector position = G4ThreeVector(0.,0.,0.);
  // point in volume test
  EInside inside = vol_1->Inside(position);
  return;
}

/*
 * ray fire test, distance to in
 * only for points external to a volume
 */
TEST_F(DagSolidTest,test_2) {
  DagMC* dagmc = DagMC::instance(); // create dag instance
  const char* h5mfilename = "/data/opt/DAGMC/Geant4/test/geometry/test_geom.h5m";
  dagmc->load_file(h5mfilename,0);
  dagmc->init_OBBTree();

  // new volume
  DagSolid* vol_1 = new DagSolid("vol_1",dagmc,1);

  G4ThreeVector position = G4ThreeVector(-100.,0.,0.);
  double distance = vol_1->DistanceToIn(position);
  std::cout << distance << std::endl;
  return;
}

/*
 * ray fire test, distance to out
 */
TEST_F(DagSolidTest,test_3) {
  DagMC* dagmc = DagMC::instance(); // create dag instance
  const char* h5mfilename = "/data/opt/DAGMC/Geant4/test/geometry/test_geom.h5m";
  dagmc->load_file(h5mfilename,0);
  dagmc->init_OBBTree();

  // new volume
  DagSolid* vol_1 = new DagSolid("vol_1",dagmc,1);

  G4ThreeVector position = G4ThreeVector(-100.,0.,0.);
  G4ThreeVector direction = G4ThreeVector(1.,0.,0.);
  double distance = vol_1->DistanceToOut(position,direction);
  std::cout << distance << std::endl;
  return;
}

/*
 * ray fire test, distance to out, called 
 *  when point is inside of volume
 */
TEST_F(DagSolidTest,test_4) {
  DagMC* dagmc = DagMC::instance(); // create dag instance
  const char* h5mfilename = "/data/opt/DAGMC/Geant4/test/geometry/test_geom.h5m";
  dagmc->load_file(h5mfilename,0);
  dagmc->init_OBBTree();

  // new volume
  DagSolid* vol_1 = new DagSolid("vol_1",dagmc,1);

  G4ThreeVector position = G4ThreeVector(0.,0.,0.);
  double distance = vol_1->DistanceToOut(position);
  std::cout << distance << std::endl;
  return;
}

/*
 * ray fire test, distance to in 
 * points inside of a volume
 */
TEST_F(DagSolidTest,test_5) {
  DagMC* dagmc = DagMC::instance(); // create dag instance
  const char* h5mfilename = "/data/opt/DAGMC/Geant4/test/geometry/test_geom.h5m";
  dagmc->load_file(h5mfilename,0);
  dagmc->init_OBBTree();

  // new volume
  DagSolid* vol_1 = new DagSolid("vol_1",dagmc,1);

  G4ThreeVector position = G4ThreeVector(0.,0.,0.);
  G4ThreeVector direction = G4ThreeVector(1.,0.,0.);
  double distance = vol_1->DistanceToIn(position,direction);
  std::cout << distance << std::endl;
  return;
}

/*
 * ray fire test, distance to out calc normal as well
 * points inside of a volume
 */
TEST_F(DagSolidTest,test_6) {
  DagMC* dagmc = DagMC::instance(); // create dag instance
  const char* h5mfilename = "/data/opt/DAGMC/Geant4/test/geometry/test_geom.h5m";
  dagmc->load_file(h5mfilename,0);
  dagmc->init_OBBTree();

  // new volume
  DagSolid* vol_1 = new DagSolid("vol_1",dagmc,1);

  // point inside cell looking out
  G4ThreeVector position = G4ThreeVector(0.,0.,0.);
  G4ThreeVector direction = G4ThreeVector(1.,0.,0.);
  
  G4ThreeVector normal;
  bool v_norm = false;
  double distance = vol_1->DistanceToOut(position,direction,true,&v_norm,&normal);
  std::cout << v_norm << std::endl;
  std::cout << normal << std::endl;
  std::cout << distance << std::endl;
  return;
}

/*
 * ray fire test, distance to out calc normal as well
 * points inside of a volume, point on surface
 */
TEST_F(DagSolidTest,test_7) {
  DagMC* dagmc = DagMC::instance(); // create dag instance
  const char* h5mfilename = "/data/opt/DAGMC/Geant4/test/geometry/test_geom.h5m";
  dagmc->load_file(h5mfilename,0);
  dagmc->init_OBBTree();

  // new volume
  DagSolid* vol_1 = new DagSolid("vol_1",dagmc,1);

  // point inside cell looking out
  G4ThreeVector position = G4ThreeVector(50.,0.,0.);
  G4ThreeVector direction = G4ThreeVector(1.,0.,0.);
  
  G4ThreeVector normal;
  bool v_norm = false;
  double distance = vol_1->DistanceToOut(position,direction,true,&v_norm,&normal);


  return;
}

/*
 * ray fire test, distance to out calc normal as well
 * point just outside of a volume, vnorm should be false
 */
TEST_F(DagSolidTest, test_8 ) {
  DagMC* dagmc = DagMC::instance(); // create dag instance
  const char* h5mfilename = "/data/opt/DAGMC/Geant4/test/geometry/test_geom.h5m";
  dagmc->load_file(h5mfilename,0);
  dagmc->init_OBBTree();

  // new volume
  DagSolid* vol_1 = new DagSolid("vol_1",dagmc,1);

  // point inside cell looking out
  G4ThreeVector position = G4ThreeVector(51.,0.,0.);
  G4ThreeVector direction = G4ThreeVector(1.,0.,0.);
  
  G4ThreeVector normal;
  bool v_norm = false;
  double distance = vol_1->DistanceToOut(position,direction,true,&v_norm,&normal);

  return;
}
