// FluDAG/src/UnitNumberManager.cpp

#include "UnitNumberManager.hpp"

/* Start and end logical (Fortran-style) unit numbers for record writing */
const int BAD_UNIT_NUMBER =   0;

const int UnitNumberManager::START_UNIT;
const int UnitNumberManager::END_UNIT;

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
UnitNumberManager::UnitNumberManager() 
{
    num_units_in_use = 0;
}

// For testing.  Make friendly, etc.
int UnitNumberManager::getNumUnitsInUse()
{
    return num_units_in_use;
}
//---------------------------------------------------------------------------//
// getUnitNumber()
//---------------------------------------------------------------------------//
// Find the unit number in a map, or calculate it and add a new map entry
int UnitNumberManager::getUnitNumber(std::string name)
{
    // During preprocessing of the string it may be determined to be invalid and not set
    if ( name.empty())
    {
       std::cout << "Warning:  Invalid name." << std::endl;
       return -1;
    }
    // If name is already in the map, return the int associated with it
    if ( UnitNumberManager::nameNumberMap.count(name)>0 )
    {
       return nameNumberMap.find(name)->second;
    }
    else // otherwise get the next unit number
    {
       int unitNumber = getNextUnitNumber();
       if (unitNumber == BAD_UNIT_NUMBER)
       {
          std::cerr << "Warning:  There are no more available unit numbers.  A unit number of 0 is being assigned." 
                    << std::endl;
       }
       else  // If the unit number is bad, no insertion is done
       {
          nameNumberMap.insert(std::make_pair(name,unitNumber)); 
       }

       return unitNumber;
    }
}

//---------------------------------------------------------------------------//
// getNextUnitNumber()
//---------------------------------------------------------------------------//
// Convenience method to get the next logical unit number for the writing-out 
// field of a FLUKA card.  The key to this function is when to call it.  
int UnitNumberManager::getNextUnitNumber()
{
    // Note that the unit number is increasingly negative, 
    // counting down from the (negative) START_UNIT
    int retval =  START_UNIT - num_units_in_use;
  
    if (retval < END_UNIT)
    {
       retval = BAD_UNIT_NUMBER;      
    }
    else  // only increase num_units_in_use if we're adding to the list
    {  
       ++num_units_in_use;
    }
   
    return retval;
}
// end of FluDAG/src/UnitNumberManager.cpp
