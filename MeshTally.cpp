#include "MeshTally.hpp" 

#include "moab/Interface.hpp"

#include <sstream>
#include <iostream>

moab::ErrorCode MeshTally::setup_tags( moab::Interface* mbi, const char* prefix ){
  
  moab::ErrorCode rval;

  std::string pfx = prefix;

  tally_tags.resize( ebins );
  error_tags.resize( ebins );

  for( unsigned i = 0; i < ebins; ++i ){

    std::string t_name = pfx + "TALLY_TAG", e_name = pfx + "ERROR_TAG";
    std::stringstream str;  

    if( i + 1 != ebins ){
      str << "_" << fmesh.energy_bin_bounds[i] << '-' << fmesh.energy_bin_bounds[i+1];
    }  
    t_name += str.str();
    e_name += str.str();

    rval = mbi->tag_get_handle( t_name.c_str(), 1, moab::MB_TYPE_DOUBLE, 
                                tally_tags[i], moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );
    if( rval != moab::MB_SUCCESS ) return rval;
    
    rval = mbi->tag_get_handle( e_name.c_str(), 1, moab::MB_TYPE_DOUBLE,
                                error_tags[i], moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );
    if( rval != moab::MB_SUCCESS ) return rval;
  }

  return moab::MB_SUCCESS;

}
