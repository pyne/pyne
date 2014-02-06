// FluDAG/src/test/test_FlukaFuncs.cpp

#include <gtest/gtest.h>

// #include "moab/Core.hpp"        // moab::Interface
#include "DagMC.hpp"
#include "MBInterface.hpp"
#include "../DagWrappers.hh"


#include <cassert>


#define DAG moab::DagMC::instance()

const std::string directory = "/filespace/people/z/zachman/repos/fludag_testing/";

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

       rval = DAG->load_file(infile.c_str(), 0.0 ); 
       assert(rval == MB_SUCCESS);

       // DAG call to initialize geometry
       rval = DAG->init_OBBTree();
       assert (rval == MB_SUCCESS);
       std::cout << "FluDAGTest:Setup():  Past init_OBBTree...." << std::endl;

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
// FIXTURE-BASED TESTS: WrapperTest
//---------------------------------------------------------------------------//
TEST_F(FluDAGTest, WrapperTest)
{
  std::cout << "Calling g_fire. " << std::endl;  
  g_fire(oldReg, point, dir, propStep, retStep, newReg);
  EXPECT_EQ(0.0, retStep);
 
}
//---------------------------------------------------------------------------//
// end of FluDAG/src/test/test_FlukaFuncs.cpp
