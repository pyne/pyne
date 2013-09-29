class TallyData
{
   public: 
   // >>> TALLY DATA ARRAYS

    /// Data array for storing sum of scores for all particle histories
    std::vector<double> tally_data;

    /// Data array for determining error in tally results
    std::vector<double> error_data;

    /// Data array for storing sum of scores for a single history
    std::vector<double> temp_tally_data;

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
     * \brief Resize data arrays to hold all of the mesh tally data
     * \param num_tally_points number of tally points included in mesh tally
     *
     * Arrays will be resized to the given number of tally points multiplied
     * by the number of energy bins.
     */
    void resize_data_arrays(unsigned int num_tally_points);
//---------------------------------------------------------------------------//
    /**
     * \brief Determines location of element in data array
     * \param data array containing element to be accessed
     * \param tally_point entity handle representing tally point
     * \param energy_bin index representing energy bin
     * \return reference to element in data array
     *
     * Enables direct access to the mesh tally data for the given tally point
     * and energy bin.
     */
    double& get_data(std::vector<double>& data,
                     moab::EntityHandle tally_point,
                     unsigned energy_bin = 0);

}
