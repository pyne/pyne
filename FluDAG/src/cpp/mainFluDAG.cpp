//----------------------------------*-C++, Fortran-*----------------------------------//
/*!
 * \file   ~/DAGMC/FluDAG/src/cpp/mainFluDAG.cpp
 * \author Julie Zachman 
 * \date   Apr 5 2013 
 * \brief  Functions called by fluka
 * \note   Unittested
 */
//---------------------------------------------------------------------------//
#include "fludag_utils.h"
#include "fluka_funcs.h"

#include "DagMC.hpp"
#include "MBInterface.hpp"
#include "MBCartVect.hpp"
#include "MBTypes.h"

#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <string.h>
#include <fstream>

using namespace moab;

#define flukam flukam_
#define DAG DagMC::instance()

extern "C" 
{ 
  void flukam(const int &GeoFlag);
}

int main(int argc, char* argv[]) 
{
  bool flukarun = false;
  MBErrorCode error;
  
  // Default h5m filename is for fluka runs
  std::string infile = "dagmc.h5m"; // used to be test.h5m should be dagmc.h5m
 
  if ( argc == 1 ) // then its a fluka run
    {
      infile = "../"+infile; // fluka create a run directory one higher than where we launch
      flukarun = true;
    }
  else if ( argc > 2 )
    {
      std::cout << "run as main_fludag <facet_file>  to produce"
	        << " material assignments" << std::endl;
      std::cout << "too many arguments provided" << std::endl;
      exit(1);
    }
  else // its a pre process run
    {
      infile = argv[1]; // must be the 2nd argument
    }


  std::ifstream h5mfile (infile.c_str()); // filestream for mesh geom
  if ( h5mfile.good() )
    {
      // ok
    }
  else // file not ok
    {
      std::cout << "h5m file does not exist" << std::endl;
      exit(1);
    }

  int max_pbl = 1;

  // load the dag file
  std::cout << "Loading the faceted geometry file " << infile << "..." << std::endl;
  error = DAG->load_file(infile.c_str(), 0.0 ); // load the dag file takeing the faceting from h5m
  if ( error != MB_SUCCESS ) 
    {
      std::cerr << "DAGMC failed to read input file: " << infile << std::endl;
      exit(EXIT_FAILURE);
    }
 
  // initialize geometry
  error = DAG->init_OBBTree();
  if ( error != MB_SUCCESS ) 
    {
      std::cerr << "DAGMC failed to initialize geometry and create OBB tree" <<  std::endl;
      exit(EXIT_FAILURE);
    }

  // if fluka preprocess run then create mat file to paste into input deck
  if (!flukarun)
    {
      std::string lcad = "mat.inp";
      fludagwrite_assignma(lcad);
      std::cout << "Producing material snippets" << std::endl;
      std::cout << "please paste these into your input deck" << std::endl;
    }
  else // call fluka run
    {
      // flugg mode is flag = 1
       const int flag = 1;
       flukam(flag);
    }

  return 0;
}

