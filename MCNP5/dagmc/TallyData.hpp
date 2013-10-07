// MCNP5/dagmc/TallyData.hpp

#ifndef DAGMC_TALLY_DATA_HPP
#define DAGMC_TALLY_DATA_HPP

#include <vector>
#include <cassert>
#include <set>
#include <utility>

class TallyData
{
   public: 
    /**
     * \brief Constructor
     */
    TallyData(unsigned int num_energy_bins, bool total_energy_bin);

    //---------------------------------------------------------------------------//
    // >>> TALLY DATA ARRAYS
    /**
     * \brief get_tally_data(), get_error_data(), get_scratch_data()
     * \param length output parameter containing size of data array
     * \return pointer to the data array
     *
     * These three methods can be useful for implementing MeshTally
     * functionality in parallel.
     */
    double* get_tally_data(int& length);
    double* get_error_data(int& length);
    double* get_scratch_data(int& length);

    /**
     * \brief Resets all of the mesh tally data arrays to zero
     */
    void zero_tally_data();
    /**
     * \brief Determines location of element in data array
     * \param tally_point index representing tally point
     * \param energy_bin index representing energy bin
     * \return reference to element in data array
     *
     * Enables direct access to the mesh tally data for the given tally point
     * and energy bin.
     */
     std::pair <double,double> get_data(unsigned int tally_point_index, unsigned int energy_bin);

     unsigned int get_num_energy_bins();
     bool has_total_energy_bin();
     void end_history();
     void add_score_to_tally(unsigned int tally_point_index, double score, unsigned int ebin);

    /**
     * \brief Resize data arrays to hold all of the mesh tally data
     * \param num_tally_points number of tally points included in mesh tally
     *
     * Arrays will be resized to the given number of tally points multiplied
     * by the number of energy bins.
     */
    void resize_data_arrays(unsigned int num_tally_points);

   private: 

    /// tally points updated in current history; cleared by end_history()
    std::set<unsigned int> visited_this_history;

    /// Data array for storing sum of scores for all particle histories
    std::vector<double> tally_data;

    /// Data array for determining error in tally results
    std::vector<double> error_data;

    /// Data array for storing sum of scores for a single history
    std::vector<double> temp_tally_data;

    /// Number of energy bins implemented in the data arrays
    unsigned int num_energy_bins;

    /// Set to true by default; determines if total energy bin is included
    bool total_energy_bin;

};

#endif // DAGMC_TALLY_DATA_HPP

// end of MCNP5/dagmc/TallyData.hpp


