// MCNP5/dagmc/TallyData.hpp

#ifndef DAGMC_TALLY_DATA_HPP
#define DAGMC_TALLY_DATA_HPP

#include <vector>
#include <set>
#include <utility>

/**
 * \class TallyData
 * \brief Defines structure for storing and accessing all tally data
 *
 * TallyData is a class that is used to store all of the Tally data that
 * is needed for use in Monte Carlo particle codes.  Three separate data
 * arrays are included with each TallyData object
 *
 *    1) tally_data: stores sum of scores for all particle histories
 *    2) error_data: stores data needed to determine error in tally results
 *    3) temp_tally_data: stores sum of scores for a single history
 *
 * Each element in these data arrays represents one tally point and one energy
 * bin.  They are ordered first by tally point, and then by energy bin.
 *
 * =======================
 * TallyData Functionality
 * =======================
 *
 * The first thing that a Tally needs to do in order to use a TallyData object
 * for storing its data is to use resize_data_arrays() based on the number of
 * tally points that are required.  The size of these arrays will then be set
 * to the number of tally points multiplied by the number of energy bins that
 * were requested when the TallyData object was constructed.
 *
 * As scores are accumulated, the Tally should then use the add_score_to_tally()
 * method for adding those scores to the appropriate data arrays.  Once a
 * particle history is complete, then the end_history() method can be used to
 * update the tally and error data arrays.
 *
 * To read tally and error values for a single tally point, the get_data()
 * function can be used.  If direct access to the underlying data structures
 * are needed, then get_tally_data(), get_error_data() and get_scratch_data()
 * can be used instead.  However, most functionality can be implemented through
 * use of other TallyData methods and direct access is not typically needed.
 */
class TallyData
{
   public: 
    /**
     * \brief Constructor
     * \param[in] num_energy_bins number of energy groups to store
     * \param[in] total_energy_bin if true, include a total energy bin
     */
    TallyData(unsigned int num_energy_bins, bool total_energy_bin);

    // >>> PUBLIC INTERFACE

    /**
     * \brief Gets tally and error values for a single tally point
     * \param[in] tally_point_index the index representing the tally point
     * \param[in] energy_bin the index representing the energy bin
     * \return pair containing (tally, error) data for the tally point
     */
     std::pair <double,double> get_data(unsigned int tally_point_index,
                                        unsigned int energy_bin) const;

    /**
     * \brief get_tally_data(), get_error_data(), get_scratch_data()
     * \param[out] length the size of the data array
     * \return pointer to the data array
     *
     * Provides direct access to the data arrays for all tally points.
     */
    double* get_tally_data(int& length);
    double* get_error_data(int& length);
    double* get_scratch_data(int& length);

    /**
     * \brief Resets all data arrays for this TallyData
     */
    void zero_tally_data();

    /**
     * \brief Resize arrays stored within this TallyData
     * \param[in] num_tally_points number of tally points to store
     *
     * Arrays will be resized to the given number of tally points multiplied
     * by the number of energy bins.
     */
    void resize_data_arrays(unsigned int num_tally_points);

    /**
     * \brief get_num_energy_bins()
     * \return Number of energy bins stored in this TallyData
     *
     * This includes the total energy bin.
     */
    unsigned int get_num_energy_bins() const;

    /**
     * \brief has_total_energy_bin()
     * \return true if a total energy bin is included in this TallyData
     */
    bool has_total_energy_bin() const;

    // >>> TALLY ACTION METHODS

    /**
     * \brief Process TallyData when a particle history is completed
     */
    void end_history();

    /**
     * \brief Add a score to this TallyData for the given tally point
     * \param[in] tally_point_index the index representing the tally point
     * \param[in] score the score to be added
     * \param[in] ebin the energy bin to which score will be added
     */
    void add_score_to_tally(unsigned int tally_point_index, double score, unsigned int ebin);

   private: 
    // Data array for storing sum of scores for all particle histories
    std::vector<double> tally_data;

    // Data array for determining error in tally results
    std::vector<double> error_data;

    // Data array for storing sum of scores for a single history
    std::vector<double> temp_tally_data;

    // tally points updated in current history; cleared by end_history()
    std::set<unsigned int> visited_this_history;

    // Number of energy bins implemented in the data arrays
    unsigned int num_energy_bins;

    // Set to true by default; determines if total energy bin is included
    bool total_energy_bin;

    // Number of tally points = tally_data.size()/num_energy_bins
    unsigned int num_tally_points;
};

#endif // DAGMC_TALLY_DATA_HPP

// end of MCNP5/dagmc/TallyData.hpp
