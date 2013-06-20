// MCNP5/dagmc/Tally.hpp

#ifndef DAGMC_TALLY_HPP
#define DAGMC_TALLY_HPP


//===========================================================================//
/**
 * \class Tally
 * \brief Defines a single tally event
 *
 * Tally is a pure virtual class that is the Observer class of the Observer pattern.
 * It is part of a Monte Carlo particle transport simulation.  
 *
 */
//===========================================================================//
class Tally
{
  public:

    /**
     * \brief Defines update
     *
     *     1) blah blah blah
     */
    virtual void update() = 0;
 
    virtual void end_history() = 0;

    virtual void write_data() = 0;
};

#endif // DAGMC_TALLY_HPP

// end of MCNP5/dagmc/Tally.hpp
