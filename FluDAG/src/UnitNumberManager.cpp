//---------------------------------------------------------------------------//
// getUnitNumber()
//---------------------------------------------------------------------------//
// Find the unit number in a map, or calculate it and add a new map entry
int UnitNumberManager::getUnitNumber(std::string name)
{
    // If name is already in the map, grab the int associated with it
    if ( UnitNumberManager::nameNumberMap.count(name)>0 )
    {
       return nameNumberMap.find(name)->second;
    }
    else // get the next unit number
    {
       int unitNumber = getNextAvailableUnitNumber();
       nameNumberMap.insert(std::make_pair(name,unitNumber)); 
       if (unitNumber == BAD_UNIT_NUMBER)
       {
          std::cerr << "Warning:  There are no more available unit numbers.  A unit number of 0 is being assigned." 
                    << std::endl;
       }
       return unitNumber;
    }
}

//---------------------------------------------------------------------------//
// getNextUnitNumber()
//---------------------------------------------------------------------------//
// Convenience method to get the next logical unit number for the writing-out 
// field of a FLUKA card.  The key is when to call.  
int UnitNumberManager::getNextUnitNumber()
{
    // Note that the unit number id increasingly negative, counting down from the (negative) START_UNIT
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
