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

  virtual void SetUp()
  {
    DagMC* dagmc = DagMC::instance(); // create dag instance
    
    // dag_volumes
    const char* h5mfilename = "test_geom.h5m";
    dagmc->load_file(h5mfilename,0);
    dagmc->init_OBBTree();
    
    // new volume
    vol_1 = new DagSolid("vol_1",dagmc,1);
  }

  protected:
  
  DagSolid* vol_1;

};

/*
 * Empty common setup function
 */
TEST_F(DagSolidTest,SetUp) {

}

/*
 * Point in volume test
 */
TEST_F(DagSolidTest,point_in_volume) {
  // sample position
  G4ThreeVector position = G4ThreeVector(0.,0.,0.);
  // point in volume test
  EInside inside = vol_1->Inside(position);
 
  EXPECT_EQ(kInside,inside);
 
  return;
}


/*
 * Point on surface
 */
TEST_F(DagSolidTest,point_on_surface) {

  // sample position
  G4ThreeVector position = G4ThreeVector(50.,0.,0.);
  // point in volume test
  EInside inside = vol_1->Inside(position);

  EXPECT_EQ(kSurface,inside);

  return;
}

/*
 * Point on surface tolerance, kCarTolerance/2.0 = 5.0e-9 
 * of a surface should return outside
 */
TEST_F(DagSolidTest,point_out_outside_tolerance) {

  // sample position
  G4ThreeVector position = G4ThreeVector(50.+(1.0e-8/2.0),0.,0.);
  // point in volume test
  EInside inside = vol_1->Inside(position);

  EXPECT_EQ(kOutside,inside);

  return;
}

/*
 * Point on surface tolerance, kCarTolerance/2.0 = 5.0e-9 
 * of a surface should return outside
 */
TEST_F(DagSolidTest,point_out_surface_tolerance) {

  // sample position
  G4ThreeVector position = G4ThreeVector(50.+(9.99e-9/2.0),0.,0.);
  // point in volume test
  EInside inside = vol_1->Inside(position);

  EXPECT_EQ(kSurface,inside);

  return;
}



/*
 * Point on surface tolerance, within kCarTolerance/2.0 = 5.0e-9 
 * of a surface should return on Surface
 */
TEST_F(DagSolidTest,point_on_surface_tolerance) {

  // sample position
  G4ThreeVector position = G4ThreeVector(50.-(9.999e-9/2.0),0.,0.);
  // point in volume test
  EInside inside = vol_1->Inside(position);

  EXPECT_EQ(kSurface,inside);

  return;
}


/*
 * Point on out of volume
 */
TEST_F(DagSolidTest,point_outside_volume) {

  // sample position
  G4ThreeVector position = G4ThreeVector(51.,0.,0.);
  // point in volume test
  EInside inside = vol_1->Inside(position);

  EXPECT_EQ(kOutside,inside);

  return;
}



/*
 * ray fire test, distance to in
 * only for points external to a volume
 */
TEST_F(DagSolidTest,test_2) {

  G4ThreeVector position = G4ThreeVector(-100.,0.,0.);
  double distance = vol_1->DistanceToIn(position);

  // DistanceToIn point external to volume should give distance
  // to entry, distance should be set to 50.0 mm 
  EXPECT_EQ(50.0,distance);

  return;
}

/*
 * ray fire test, distance to out
 */
TEST_F(DagSolidTest,test_3) {

  G4ThreeVector position = G4ThreeVector(-100.,0.,0.);
  G4ThreeVector direction = G4ThreeVector(1.,0.,0.);
  double distance = vol_1->DistanceToOut(position,direction);

  // DistanceToOut normally only called from inside the volume
  // but in this case, we expect the ray to skip over the volume 
  // and give the distance to when we leave the volume instead

  // distance should be set to 50.0 mm since the point
  // is in the center of the box
  EXPECT_EQ(150.0,distance);

  return;
}

/*
 * ray fire test, distance to out, called 
 *  when point is inside of volume
 */
TEST_F(DagSolidTest,test_4) {

  G4ThreeVector position = G4ThreeVector(0.,0.,0.);
  double distance = vol_1->DistanceToOut(position);

  // distance should be set to 50.0 mm since the point
  // is in the center of the box
  EXPECT_EQ(50.0,distance);

  return;
}

/*
 * ray fire test, distance to in 
 * points inside of a volume
 */
TEST_F(DagSolidTest,test_5) {

  G4ThreeVector position = G4ThreeVector(0.,0.,0.);
  G4ThreeVector direction = G4ThreeVector(1.,0.,0.);
  double distance = vol_1->DistanceToIn(position,direction);

  // distance should be set to infinity since there are no other volumes
  // to enter
  EXPECT_EQ(kInfinity,distance);

  return;
}

/*
 * ray fire test, distance to out calc normal as well
 * points inside of a volume
 */
TEST_F(DagSolidTest,test_6) {

  // point inside cell looking out
  G4ThreeVector position = G4ThreeVector(0.,0.,0.);
  G4ThreeVector direction = G4ThreeVector(1.,0.,0.);
  
  G4ThreeVector normal;
  bool v_norm = false;
  double distance = vol_1->DistanceToOut(position,direction,true,&v_norm,&normal);

  // when the point is inside or on the surface of the volume, v_norm should be true
  EXPECT_TRUE(v_norm);
  // distance should be set to 0.0 in this case
  EXPECT_EQ(50.0,distance);
  
  return;
}

/*
 * ray fire test, distance to out calc normal as well
 * points inside of a volume, point on surface
 */
TEST_F(DagSolidTest,test_7) {

  // point inside cell looking out
  G4ThreeVector position = G4ThreeVector(50.,0.,0.);
  G4ThreeVector direction = G4ThreeVector(1.,0.,0.);
  
  G4ThreeVector normal;
  bool v_norm = false;
  double distance = vol_1->DistanceToOut(position,direction,true,&v_norm,&normal);

  // when the point is inside or on the surface of the volume, v_norm should be true
  EXPECT_TRUE(v_norm);
  // distance should be set to 0.0 in this case
  EXPECT_EQ(0.0,distance);
  return;
}

/*
 * ray fire test, distance to out calc normal as well
 * point just outside of a volume, vnorm should be false
 */
TEST_F(DagSolidTest, test_8 ) {

  // point inside cell looking out
  G4ThreeVector position = G4ThreeVector(51.,0.,0.);
  G4ThreeVector direction = G4ThreeVector(1.,0.,0.);
  
  G4ThreeVector normal;
  bool v_norm = false;
  double distance = vol_1->DistanceToOut(position,direction,true,&v_norm,&normal);
 
  // when the point is outside the volume, v_norm should be false
  EXPECT_FALSE(v_norm);
  // distance should be set to infinity
  EXPECT_EQ(kInfinity,distance);
  return;
}


/*
 * volume_test calculates the volume of the solid, test cube is 10*10*10 cm
 * G4 works in mm, therefore expect 10*10*10*1000 = 1e6 cubic millimetres
 */
TEST_F(DagSolidTest, volume_test ) {

  G4double volume =  vol_1->GetCubicVolume();
  EXPECT_EQ(1.0e6,volume);
  return;
}


/*
 * surface_area_test determines if DAGMC returns correct surface area
 * surface area currently is set to zero
 */
TEST_F(DagSolidTest, surface_area_test ) {

  G4double surface_area = vol_1->GetSurfaceArea();
  std::cout << surface_area << std::endl;
  return;
}


