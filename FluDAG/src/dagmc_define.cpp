// FluDAG/src/dagmc_define.cpp
//----------------------------------*-C++, Fortran-*----------------------------------//
/*!
 * \file   ~/DAGMC/FluDAG/src/dagmc_define.cpp
 * \author Julie Zachman 
 * \date   Dec 10 2013 
 * \brief  Functions called by fluka
 * \note   Unittested
 */
//---------------------------------------------------------------------------//
#include "dagmc_define.hpp"  // fluka_funcs
#include "chkerr.hpp"

#include "moab/Interface.hpp"
#include "moab/Core.hpp"
#include "moab/ProgOptions.hpp"

#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <string.h>
#include <fstream>
#include <time.h>       // for timing the routine

using namespace moab;

  // void fludagwrite_assignma(moab::Interface* mbi, std::string matfile)
   //{
    //   return;  
   //}

int main(int argc, char* argv[]) 
{
  ProgOptions po("dagmc_define: a tool for processing h5m files for physics analysis");
  std::string dagversion;
  // DagMC::version( &dagversion );
  // po.setVersion( dagversion );

  std::string input_file;
  std::string output_file = "dagmc_define_out.h5m";
  bool verbose = false;
  po.addOpt<void>( ",v", "Verbose output", &verbose );
  // po.addOpt<std::string>( "outmesh,o", "Specify output file name (default "+output_file+")", &output_file );
  // po.addOpt<void>("no-outmesh,", "Do not write an output mesh" );
  // po.addOpt<std::string>( ",m", "Specify alternate input mesh to override surfaces in input_file" );
  // po.addOpt<std::string>( "obb-vis,O", "Specify obb visualization output file (default none)" );
  // po.addOpt<int>( "obb-vis-divs", "Resolution of obb visualization grid (default 50)", &grid );
  // po.addOpt<void>( "obb-stats,S", "Print obb statistics.  With -v, print verbose statistics." );
  // po.addOpt<std::vector<int> >( "vols,V", "Specify a set of volumes (applies to --obb_vis and --obb_stats, default all)" );
  po.addOpt<std::string>("fluka,", "Specify create FlUKA records.");
  po.addOpt<std::string>("geant,", "Specify create GEANT4 records and code.");
  po.addOpt<std::string>("input,i", "Specify the .h5m file for processing", &input_file);
  po.addOpt<std::string>("output,o", "Specify the output file", &output_file);

  time_t time_before,time_after;

  // Default h5m filename is for fluka runs
  // std::string infile = "dagmc.h5m"; 
 
  moab::Interface* mbi = new moab::Core();
/*
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
*/

  std::ifstream h5mfile (input_file.c_str()); // filestream for mesh geom
  if ( !h5mfile.good() )
  {
      // file not ok
      std::cout << "input h5m file does not exist" << std::endl;
      exit(1);
   }

  // int max_pbl = 1;

  // get the current time
  time(&time_before);  /* get current time; same as: timer = time(NULL)  */


  // load the dag file
  std::cout << "Loading the faceted geometry file " << input_file << "..." << std::endl;

  // error = DAG->load_file(infile.c_str(), 0.0 ); // load the dag file takeing the faceting from h5m
  // error = mbi
 //  if ( error != MB_SUCCESS ) 
   // {
    //  std::cerr << "DAGMC failed to read input file: " << infile << std::endl;
     // exit(EXIT_FAILURE);
    //}
  
  EntityHandle input_file_set;
  ErrorCode error;

  error = mbi->create_meshset( 0, input_file_set );
  CHECKERR( mbi, error );

  error = mbi->load_file( input_file.c_str(), &input_file_set );
  CHECKERR( mbi, error );

  // Timing report
  time(&time_after); 
  double seconds = difftime(time_after,time_before); //get the time in seconds to load file
  time_before = time_after; // reset time to now for the next call
  std::cout << "Time to load the h5m file = " << seconds << " seconds" << std::endl;

  // ToDo:  Do I need to do this?
  // initialize geometry
/*
  error = DAG->init_OBBTree();
  if ( error != MB_SUCCESS ) 
    {
      std::cerr << "DAGMC failed to initialize geometry and create OBB tree" <<  std::endl;
      exit(EXIT_FAILURE);
    }
*/
  time(&time_after);

  seconds = difftime(time_after,time_before);
  std::cout << "Time to initialise the geometry" << seconds << std::endl;

  // Main call to parse the file and write out the records
  std::string lcad = "mat.inp";
  fludagwrite_assignma(mbi, lcad);
  std::cout << "Producing material snippets" << std::endl;
  std::cout << "please paste these into your input deck" << std::endl;

  return 0;
}

