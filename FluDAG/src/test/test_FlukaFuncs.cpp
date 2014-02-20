// FluDAG/src/test/test_FlukaFuncs.cpp

#include <gtest/gtest.h>

// #include "moab/Core.hpp"        // moab::Interface
#include "DagMC.hpp"
#include "MBInterface.hpp"
#include "../DagWrappers.hh"


#include <cassert>


#define DAG moab::DagMC::instance()

const std::string directory = "/filespace/people/z/zachman/repos/fludag_testing/";
int num_slab_vols = 12;

//---------------------------------------------------------------------------//
// HELPER METHODS
//---------------------------------------------------------------------------//
// loads default mesh into the mbi instance and gets all mesh nodes
/*
void load_default_mesh(moab::Interface* mbi, moab::Range& mesh_nodes)
{
    std::string default_slabs = dir + "input/test2/fludag/slabs.h5m";

    moab::ErrorCode rval = mbi->load_mesh(default_slabs);
    assert(rval == moab::MB_SUCCESS);

    // copy all mesh nodes into Range
    moab::EntityHandle root_set = 0;
    rval = mbi->get_entities_by_type(root_set, moab::MBVERTEX, mesh_nodes);
    assert(rval == moab::MB_SUCCESS);
}
*/
//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//
class FluDAGTest : public ::testing::Test
{
  protected:

    // initalize variables for each test
    virtual void SetUp()
    {
       // Default h5m file for testing
       // std::string infile = directory + "input/test2/fludag/slabs.h5m";
       std::string infile = "../slabs.h5m";

       rloadval = DAG->load_file(infile.c_str(), 0.0 ); 
       assert(rval == MB_SUCCESS);

       // DAG call to initialize geometry
       rval = DAG->init_OBBTree();
       assert (rval == MB_SUCCESS);

       // Initialize point and dir
       point[0] = 0.0;
       point[1] = 0.0; 
       point[2] = 0.0;

       dir[0] = 0.0;
       dir[1] = 0.0; 
       dir[2] = 0.0;
  
       oldReg = 1;
       propStep = 0.0;
    }

  protected:

    MBErrorCode rloadval;
    MBErrorCode rval;

    // Position
    double point[3];

    // Direction Vector
    double dir[3];
    
    int      oldReg;
    double propStep;

    double retStep;
    int    newReg;
   
//---------------------------------------------------------------------------//
// void g_fire(int& oldRegion, double point[], double dir[], 
//              double &propStep, double& retStep,  int& newRegion)
//---------------------------------------------------------------------------//
// oldRegion - the region of the particle's current coordinates
// point     - the particle's coordinate location vector
// dir       - the direction vector of the particle's current path (ray)
// propStep  - ??
// retStep   - returned as the distance from the particle's current location, along its ray, to the next boundary
// newRegion - gotten from the value returned by DAG->next_vol
// newRegion is gotten from the volue returned by DAG->next_vol
// void g_fire(int& oldRegion, double point[], double dir[], double &propStep, double& retStep,  int& newRegion)

};
//---------------------------------------------------------------------------//
// Test setup outcomes
TEST_F(FluDAGTest, SetUp)
{
    EXPECT_EQ(MB_SUCCESS, rloadval);

    // DAG call to initialize geometry
    EXPECT_EQ(MB_SUCCESS, rval);

    int num_vols = DAG->num_entities(3);
    std::cout << "Number of regions is " << num_vols << std::endl;
    EXPECT_EQ(num_slab_vols, num_vols);

    std::vector< std::string > keywords;
    rval = DAG->detect_available_props( keywords );
    EXPECT_EQ(MB_SUCCESS, rval);
    rval = DAG->parse_properties( keywords );
    EXPECT_EQ(MB_SUCCESS, rval);
    
    int ret, volume;
    for (unsigned i=1; i<=num_slab_vols; i++)
    {
      int id = DAG->id_by_index(3, i);
      std::cout << "Vol " << i << ", id = " << id << std::endl; 

      moab::EntityHandle eh = DAG->entity_by_index(3,i);
      rval = DAG->point_in_volume(eh, point, ret); 
      EXPECT_EQ(MB_SUCCESS, rval);
      if (ret == 1)
      {
         volume = i;
         std::cout << "\tPoint is in this volume!" << std::endl;
      } 
    }
}
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: WrapperTest
//---------------------------------------------------------------------------//
TEST_F(FluDAGTest, GFireTests)
{
  std::cout << "Calling g_fire. Start in middle of leftmost cube" << std::endl;  
  oldReg   = 2;
  point[2] = 5.0;
  dir[2]   = 1.0;
  
  g_fire(oldReg, point, dir, propStep, retStep, newReg);
  EXPECT_EQ(5.0, retStep);
}
//---------------------------------------------------------------------------//
// end of FluDAG/src/test/test_FlukaFuncs.cpp
