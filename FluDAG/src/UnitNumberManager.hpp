// FluDAG/src/UnitNumberManager.hpp

#ifndef FLUDAG_SRC_UNIT_NUMBER_MANAGER_HPP
#define FLUDAG_SRC_UNIT_NUMBER_MANAGER_HPP

#include <iostream>   // std::cerr
#include <string>     // std::string
#include <map>        // std::map

// NOTE:  the logical unit numbers in the Fluka cards are floats
class UnitNumberManager
{
  public:
    /**
     * \brief Constructor
     */
    UnitNumberManager();


    int getNumUnitsInUse();
    int getUnitNumber(std::string name);

    // Need access for test verification
    static const int START_UNIT = -21;
    static const int END_UNIT   = -99;

  private:

    std::map<std::string, int> nameNumberMap;
    int getNextUnitNumber();
    /* 
     * The number of types of cards to be written. Card 'types' are distinguished by
     * the combination of particle and tally type 
     */
    int num_units_in_use;
};

#endif // FLUDAG_SRC_UNIT_NUMBER_MANAGER_HPP

// end of FluDAG/src/UnitNumberManager.hpp
