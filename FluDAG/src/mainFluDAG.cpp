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
#include <time.h>       // for timing the routine

using namespace moab;

#define flukam flukam_
#define DAG DagMC::instance()

extern "C" 
{ 
  void flukam(const int &GeoFlag);
}
// This has been modified so that all runs are now fluka runs.  
// Getting the geometry out of the file is done by dagmc_define
int main(int argc, char* argv[]) 
{
  MBErrorCode error;
  time_t time_before,time_after;

  // Default h5m filename is for fluka runs
  std::string infile = "dagmc.h5m"; 
 
  if ( argc >= 1 ) // then its a fluka run
    {
      // fluka creates a run dir one lvl higher 
      infile = "../"+infile; 
    }
  else if ( argc > 2 )
    {
      std::cout << "run as main_fludag <facet_file>  to produce"
	        << " material assignments" << std::endl;
      std::cout << "too many arguments provided" << std::endl;
      exit(1);
    }
//  else // its a pre process run
//    {
//      infile = argv[1]; // must be the 2nd argument
//    }


  std::ifstream h5mfile (infile.c_str()); // filestream for mesh geom
  if ( !h5mfile.good() )
  {
      std::cout << "h5m file does not exist" << std::endl;
      exit(1);
  }

  int max_pbl = 1;

  // get the current time
  time(&time_before);  /* get current time; same as: timer = time(NULL)  */


  // DAG call to load the file
  std::cout << "Loading the faceted geometry file " << infile << "..." << std::endl;
  error = DAG->load_file(infile.c_str(), 0.0 ); // load the dag file takeing the faceting from h5m
  if ( error != MB_SUCCESS ) 
    {
      std::cerr << "DAGMC failed to read input file: " << infile << std::endl;
      exit(EXIT_FAILURE);
    }
  
  time(&time_after); 

  double seconds = difftime(time_after,time_before); //get the time in seconds to load file

  time_before = time_after; // reset time to now for the next call
  
  std::cout << "Time to load the h5m file = " << seconds << " seconds" << std::endl;

  // DAG call to initialize geometry
  error = DAG->init_OBBTree();
  if ( error != MB_SUCCESS ) 
    {
      std::cerr << "DAGMC failed to initialize geometry and create OBB tree" <<  std::endl;
      exit(EXIT_FAILURE);
    }

  time(&time_after);

  seconds = difftime(time_after,time_before);
  std::cout << "Time to initialise the geometry" << seconds << std::endl;

  // flugg mode is flag = 1
  const int flag = 1;
  flukam(flag);

  return 0;
}

