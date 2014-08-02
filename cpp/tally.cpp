// Tally.cpp
// Central Tally Class
// -- Andrew Davis

#include <string>
#include <vector>

#ifndef PYNE_IS_AMALGAMATED
  #include "tally.h"
#endif


/***************************/
/*** Protected Functions ***/
/***************************/

/// there are no protected functions currently
/// fool.

/************************/
/*** Public Functions ***/
/************************/

/*--- Constructors ---*/
pyne::Tally::Tally()
{
  // Empty Tally Constructor
  tally_type = "";
  particle_name = "";
  entity_id = -1;
  entity_type = "";
  entity_name = "";
}

pyne::Tally::Tally(std::string type, std::string part_name, 
		   int ent, std::string ent_type, 
		   std::string ent_name)
{
  // Empty Tally Constructor
  tally_type = type;
  particle_name = part_name;
  entity_id = ent;
  entity_type = ent_type;
  entity_name = ent_name;
}

// Destructor
pyne::Tally::~Tally()
{
};

/*--- Method definitions ---*/

std::ostream& operator<<(std::ostream& os, pyne::Tally tal)
{
  //print the Tally to ostream
  os << "\t---------\n";
  os << "\t Tallying " << tal.particle_name << " " << tal.tally_type << "\n";
  os << "\t in/on " << tal.entity_type << " " << tal.entity_id << "\n";
  return os;
};
