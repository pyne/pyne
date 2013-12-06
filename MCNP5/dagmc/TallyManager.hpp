// MCNP5/dagmc/TallyManager.hpp

#ifndef DAGMC_TALLY_MANAGER_HPP
#define DAGMC_TALLY_MANAGER_HPP

#include "Tally.hpp"
#include "TallyEvent.hpp"

//===========================================================================//
/**
 * \class TallyManager
 * \brief Defines an interface for managing all DAGMC tallies
 *
 * TallyManager is the interface that allows any Monte Carlo physics code to
 * use DAGMC tallies instead of their native implementations.  DAGMC tallies
 * can be added and removed using addNewTally() and removeTally() respectively.
 * TallyManager maintains a list of all currently active tallies, each defined
 * by a unique tally ID number.
 *
 * DAGMC tally types that are currently available include
 *
 *     "unstr_track": Unstructured tracklength mesh tally (TrackLengthMeshTally)
 *     "kde_coll": KDE collision mesh tally (KDEMeshTally)
 *     "kde_subtrack": KDE sub-track mesh tally (KDEMeshTally)
 *     "kde_track": KDE integral-track mesh tally (KDEMeshTally)
 *
 * See the individual implementations for a more detailed description and a
 * list of all available tally options for that particular tally type.  The
 * tally options are all stored as key-value pairs in a multimap, which may
 * be empty if no options are requested.
 *
 * NOTE: The DAGMC Tally implementation is based on the Observer pattern.
 * TallyManager acts as the Subject/Observable, whereas Tally objects act as
 * the Observers.  All Tally actions are performed through the TallyManager.
 *
 * ==========================
 * TallyManager Functionality
 * ==========================
 *
 * Once a list of DAGMC tallies has been created, there are a series of actions
 * that can then be performed to compute the scores and write the results for
 * all currently active tallies.  The first step is to set an event type using
 * either set_collision_event() or set_track_event().  This can be done during
 * the Monte Carlo simulation as each event occurs, or during post-processing
 * if the complete particle history data is available.
 *
 * After an event type has been set, the TallyManager can then be used to
 * updateTallies().  This will compute the scores for all currently active
 * tallies, based on the tally event data that was set.  If a particular event
 * type is not compatible with a specific DAGMC tally type, then nothing will
 * happen and the code will move on to the next tally in the list.  Note that
 * when updateTallies() has updated all of the tallies it will then reset the
 * event data using clear_last_event().
 *
 * As each particle history is completed, the end_history() method should be
 * called through the TallyManager.  This adds the current sum of scores to
 * the total for each tally that is currently active.  It is also important
 * for calculating the relative standard errors in the final tally results,
 * and may perform other functionality as needed based on the tally type.
 *
 * When all particle histories have been scored, the final step is to use the
 * write_data() method to write all of the tally results and their corresponding
 * relative standard errors to an output file.  These results are typically
 * normalized by the number of histories reported by the physics code.
 *
 * =================
 * Tally Multipliers
 * =================
 *
 * TallyManager also includes optional support for energy-dependent tally
 * multipliers.  These multipliers are defined by a unique ID determined by
 * the user, and are stored separately to the list of Tally objects so that
 * multiple tallies can use the same multiplier ID.
 *
 * To use multipliers, simply add them using addNewMultiplier() and then assign
 * them to a specific Tally using addMultiplierToTally().  This will store the
 * multiplier ID in the Tally so that it has access to that multiplier during
 * the transport process for computing its scores.  As the multiplier values
 * change, use updateMultiplier() to update their values in the TallyManager.
 */
//===========================================================================//
class TallyManager
{
  public:
    /**
     * \brief Constructor
     */
    TallyManager();
    
    // >>> PUBLIC INTERFACE

    /**
     * \brief Create a new DAGMC Tally and add it to the Observer list
     * \param[in] tally_id the unique ID for this Tally
     * \param[in] tally_type the type of Tally to create
     * \param[in] energy_bin_bounds the boundaries of the energy bins
     * \param[in] options the set of options requested for this Tally
     *
     * If an invalid tally_type is requested, then it will be ignored.
     * Note that energy_bin_bounds must have at least two entries, so
     * to define a single energy bin use [0.0, max_energy].
     */
    void addNewTally(unsigned int tally_id,
                     std::string tally_type,
                     const std::vector<double>& energy_bin_bounds,
                     std::multimap<std::string, std::string>& options);

    /**
     * \brief Add a new tally multiplier
     * \param[in] multiplier_id the unique ID for the multiplier
     *
     * All new multipliers that are added are set to a default value of 1.0.
     * If the multiplier ID already exists, then nothing will happen.
     */
    void addNewMultiplier(unsigned int multiplier_id);

    /**
     * \brief Assign a multiplier to a Tally
     * \param[in] multiplier_id the unique ID for the multiplier
     * \param[in] tally_id the unique ID for the Tally
     *
     * Note that if a multiplier ID has already been set for the Tally, then it
     * will be overwritten.  Currently only single multipliers are supported.
     */
    void addMultiplierToTally(unsigned int multiplier_id, unsigned int tally_id);

    /**
     * \brief Update the value associated with the multiplier ID
     * \param[in] multiplier_id the unique ID for the multiplier
     * \param[in] value the value of the multiplier
     *
     * If the multiplier_id is invalid, then nothing will happen.
     */
    void updateMultiplier(unsigned int multiplier_id, double value);

    /**
     * \brief Remove a DAGMC Tally from the Observer list
     * \param[in] tally_id the unique ID for the Tally to be removed
     */
    void removeTally(unsigned int tally_id);

    /**
     * \brief Set a collision event
     * \param[in] x, y, z coordinates of the collision point
     * \param[in] particle_energy the energy of the particle prior to collision
     * \param[in] particle_weight the weight of the particle prior to collision
     * \param[in] total_cross_section the macroscopic cross section for current cell
     * \param[in] cell_id the unique ID for the current cell
     * \return true if a collision event was set; false otherwise
     */
    bool setCollisionEvent(double x, double y, double z,
                           double particle_energy, double particle_weight,
                           double total_cross_section, int cell_id); 

    /**
     * \brief Set a track event
     * \param[in] x, y, z coordinates of the start of the track
     * \param[in] u, v, w current direction of the particle
     * \param[in] particle_energy the energy of the particle prior to event
     * \param[in] particle_weight the weight of the particle prior to event
     * \param[in] track_length the length of the track
     * \return true if a track event was set; false otherwise
     */
    bool setTrackEvent(double x, double y, double z,
                       double u, double v, double w,                           
                       double particle_energy, double particle_weight,
                       double track_length, int cell_id); 

    /**
     *  \brief Reset a tally event
     *
     *  Sets event type to NONE and clears all event data.
     */
    void clearLastEvent();

    /**
     * \brief Call compute_score() for all active DAGMC tallies
     *
     * Resets the tally event once all scores are computed.
     */
    void updateTallies();

    /**
     * \brief Call end_history() for all active DAGMC tallies
     */
    void endHistory();

    /**
     * \brief Call write_data() for all active DAGMC tallies
     * \param[in] num_histories the number of particle histories tracked
     */
    void writeData(double num_histories);

    // >>> TALLY DATA ACCESS METHODS

    /**
     * \brief get_tally_data(), get_error_data(), get_scratch_data()
     * \param[in] tally_id the unique ID of the Tally for which data is needed
     * \param[out] length the size of the data array
     * \return pointer to the data array
     *
     * In general, the data arrays stored in each Tally should not be accessed
     * directly.  However, if direct access to these data arrays are needed for
     * parallel or other implementations, then these methods can be used to get
     * a pointer to the underlying data stored within each Tally.
     */
    double* getTallyData(int tally_id, int& length);
    double* getErrorData(int tally_id, int& length);
    double* getScratchData(int tally_id, int& length);

    /**
     * \brief Resets all data arrays for all active Tally Observers
     */
    void zeroAllTallyData();

  private:
    // Keep a record of the currently active Tally Observers
    std::map<int, Tally*> observers; 

    // Store event data read by all active DAGMC tallies
    TallyEvent event;

    // >>> PRIVATE METHODS

    /**
     * \brief Create a new DAGMC Tally
     * \param[in] tally_id the unique ID for this Tally
     * \param[in] tally_type the type of Tally to create
     * \param[in] energy_bin_bounds the boundaries of the energy bins
     * \param[in] options the set of options requested for this Tally
     *
     * Sets up TallyInput and calls the Tally factory method.
     */
    Tally *createTally(unsigned int tally_id,
                       std::string  tally_type,
                       const std::vector<double>& energy_bin_bounds,
                       std::multimap<std::string, std::string>& options);

    /**
     * \brief Sets up TallyEvent
     * \param[in] x, y, z the position of the particle
     * \param[in] u, v, w current direction of the particle
     * \param[in] particle_energy the energy of the particle prior to event
     * \param[in] particle_weight the weight of the particle prior to event
     * \param[in] track_length the length of the track
     * \param[in] total_cross_section the macroscopic cross section for current cell
     * \param[in] cell_id the unique ID for the current geometric cell
     * \return true if an event was set; false otherwise
     */
    bool setEvent(TallyEvent::EventType type,
                   double x, double y, double z,
                   double u, double v, double w,                           
                   double particle_energy, double particle_weight,
                   double track_length, double total_cross_section,
                   int cell_id); 
};

#endif // DAGMC_TALLY_MANAGER_HPP

// end of MCNP5/dagmc/TallyManager.hpp
