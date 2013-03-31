// MCNP5/dagmc/MeshTally.cpp

#include <iostream>
#include <sstream>

#include "moab/Interface.hpp"

#include "MeshTally.hpp" 

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
MeshTally::MeshTally(const MeshTallyInput& input)
    : input_data(input)
{
    // Determine the total number of energy bins requested
    num_energy_bins = input_data.energy_bin_bounds.size();

    if(!input_data.total_energy_bin)
        --num_energy_bins;

    assert(num_energy_bins > 0);
}
//---------------------------------------------------------------------------//
// TALLY DATA ACCESS FUNCTIONS
//---------------------------------------------------------------------------//
double* MeshTally::get_tally_data(int& length)
{
    length = tally_data.size();
    return &(tally_data[0]);
}
//---------------------------------------------------------------------------//
double* MeshTally::get_error_data(int& length)
{
    length = error_data.size();
    return &(error_data[0]);
}
//---------------------------------------------------------------------------//
double* MeshTally::get_scratch_data(int& length)
{
    length = temp_tally_data.size();
    return &(temp_tally_data[0]);
}
//---------------------------------------------------------------------------//
void MeshTally::zero_tally_data()
{
    std::fill(tally_data.begin(), tally_data.end(), 0);
    std::fill(error_data.begin(), error_data.end(), 0);
    std::fill(temp_tally_data.begin(), temp_tally_data.end(), 0);
}
//---------------------------------------------------------------------------//
// PROTECTED FUNCTIONS
//---------------------------------------------------------------------------//
void MeshTally::resize_data_arrays(unsigned int num_tally_points)
{
    int new_size = num_tally_points * num_energy_bins;

    tally_data.resize(new_size, 0);
    error_data.resize(new_size, 0);
    temp_tally_data.resize(new_size, 0);
}
//---------------------------------------------------------------------------//
unsigned int MeshTally::get_entity_index(moab::EntityHandle tally_point)
{
    unsigned int ret = tally_points.index(tally_point);
    assert(ret < tally_points.size());
    return ret;
}
//---------------------------------------------------------------------------//
double& MeshTally::get_data(std::vector<double>& data,
                            moab::EntityHandle tally_point,
                            unsigned energy_bin)
{
    assert(energy_bin < num_energy_bins);
    int index = get_entity_index(tally_point) * num_energy_bins + energy_bin;
    return data[index];
}
//---------------------------------------------------------------------------//
moab::ErrorCode MeshTally::setup_tags(moab::Interface* mbi,
                                      const char* prefix)
{ 
    moab::ErrorCode rval;
    std::string pfx = prefix;
    int tag_size = 1;

    tally_tags.resize(num_energy_bins);
    error_tags.resize(num_energy_bins);

    // Create separate MOAB tag handles for every energy bin
    for(unsigned i = 0; i < num_energy_bins; ++i)
    {
        std::string t_name = pfx + "TALLY_TAG", e_name = pfx + "ERROR_TAG";
        std::stringstream str;  

        if(i + 1 != num_energy_bins)
        {
            str << "_" << input_data.energy_bin_bounds[i] 
                << '-' << input_data.energy_bin_bounds[i+1];
        }
  
        t_name += str.str();
        e_name += str.str();

        rval = mbi->tag_get_handle(t_name.c_str(),
                                   tag_size,
                                   moab::MB_TYPE_DOUBLE, 
                                   tally_tags[i],
                                   moab::MB_TAG_DENSE|moab::MB_TAG_CREAT);

        if(rval != moab::MB_SUCCESS) return rval;
    
        rval = mbi->tag_get_handle(e_name.c_str(),
                                   tag_size,
                                   moab::MB_TYPE_DOUBLE,
                                   error_tags[i],
                                   moab::MB_TAG_DENSE|moab::MB_TAG_CREAT);

        if(rval != moab::MB_SUCCESS) return rval;
    }

    return moab::MB_SUCCESS;
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/MeshTally.hpp
