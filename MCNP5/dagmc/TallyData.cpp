// MCNP5/dagmc/TallyData.cpp

#include <cassert>

#include "TallyData.hpp"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
TallyData::TallyData(unsigned int num_energy_bins, bool total_energy_bin) 
{
    assert(num_energy_bins >= 1);
    if (num_energy_bins == 1)
    {
       // Any time there is only one energy bin, total_energy_bin is moot
       this->total_energy_bin = false;
    }
    else
    {
       // Otherwise, pass through the value
       this->total_energy_bin = total_energy_bin;
    }

    //////////////////////////////////////////////////
    if (this->total_energy_bin)
    {
       // Add an extra bin for the total, if required
       this->num_energy_bins = num_energy_bins + 1;
    }
    else
    {
       this->num_energy_bins = num_energy_bins;
    }
    
    this->num_tally_points = 0;
}
//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
std::pair <double,double> TallyData::get_data(unsigned int tally_point_index,
                                              unsigned int energy_bin) const
{
   assert(energy_bin < num_energy_bins);
   assert(tally_point_index < num_tally_points);

   int index = tally_point_index * num_energy_bins + energy_bin;
   double tally = tally_data.at(index);
   double error = error_data.at(index);
 
   return std::make_pair(tally, error);
}
//---------------------------------------------------------------------------//
double* TallyData::get_tally_data(int& length)
{
    assert(tally_data.size() != 0);
    length = tally_data.size();
    return &(tally_data[0]);
}
//---------------------------------------------------------------------------//
double* TallyData::get_error_data(int& length)
{
    assert(error_data.size() != 0);
    length = error_data.size();
    return &(error_data[0]);
}
//---------------------------------------------------------------------------//
double* TallyData::get_scratch_data(int& length)
{
    assert(temp_tally_data.size() != 0);
    length = temp_tally_data.size();
    return &(temp_tally_data[0]);
}
//---------------------------------------------------------------------------//
void TallyData::zero_tally_data()
{
    std::fill(tally_data.begin(), tally_data.end(), 0);
    std::fill(error_data.begin(), error_data.end(), 0);
    std::fill(temp_tally_data.begin(), temp_tally_data.end(), 0);
}
//---------------------------------------------------------------------------//
void TallyData::resize_data_arrays(unsigned int tally_points)
{
    assert(tally_points > 0);
    num_tally_points = tally_points;
    unsigned int new_size = num_tally_points * num_energy_bins;

    tally_data.resize(new_size, 0);
    error_data.resize(new_size, 0);
    temp_tally_data.resize(new_size, 0);
}
//---------------------------------------------------------------------------//
unsigned int TallyData::get_num_energy_bins() const
{
    return num_energy_bins;
}
//---------------------------------------------------------------------------//
bool TallyData::has_total_energy_bin() const
{
    return total_energy_bin;
}
//---------------------------------------------------------------------------//
// TALLY ACTION METHODS
//---------------------------------------------------------------------------//
void TallyData::end_history()
{
    std::set<unsigned int>::iterator it;

    // add sum of scores for this history to mesh tally for each tally point
    for (it = visited_this_history.begin(); it != visited_this_history.end(); ++it)
    {
        for (unsigned int j = 0; j < num_energy_bins; ++j)
        {
            int index = (*it) * num_energy_bins + j;
            double& history_score = temp_tally_data.at(index);
            double& tally         = tally_data.at(index); 
            double& error         = error_data.at(index);

            tally += history_score;
            error += history_score * history_score;

            // reset temp_tally_data array for the next particle history
            history_score = 0;
        }
    }

    // reset set of tally points for next particle history
    visited_this_history.clear();
}
//---------------------------------------------------------------------------//
void TallyData::add_score_to_tally(unsigned int tally_point_index,
                                   double score,
                                   unsigned int energy_bin)
{
    assert(tally_point_index < num_tally_points);
    assert(energy_bin < num_energy_bins);

    // update tally for this history with new score
    int index = tally_point_index * num_energy_bins + energy_bin;;
    temp_tally_data.at(index) += score; 

    // also update total energy bin tally for this history if one exists
    if (total_energy_bin)
    {
        index = tally_point_index * num_energy_bins + num_energy_bins - 1;
        temp_tally_data.at(index) += score;
    }

    visited_this_history.insert(tally_point_index);
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/TallyData.cpp
