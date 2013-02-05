#include "DagMC.hpp"
#include <iostream>
#include <stdlib.h>
#include <string>
#include "moab/Types.hpp"

using namespace moab;

#define DAG DagMC::instance()

ErrorCode readVol(char *fileptr)
{
  // See if this works
  int num_vols = DAG->num_entities(3);
  std::cout << __FILE__ << ", " << __func__ << ":" << __LINE__ << "_______________" << std::endl;
  std::cout << "\tnum_vols in " << fileptr << " is " << num_vols << std::endl;

  ErrorCode code;
  double volume_measure;
  
  // Kerry's way:  use iterator on Range vols to get EntityHandle vols
  DagMC& dagmc = *DagMC::instance();
  Interface& moab = *dagmc.moab_instance();

  Tag dim_tag = dagmc.geom_tag();
  Range vols;
  const int three = 3;
  const void* ptr = &three;
  code = moab.get_entities_by_type_and_tag( 0, MBENTITYSET, &dim_tag, &ptr, 1, vols);
  std::cout << __FILE__ << ", " << __func__ << ":" << __LINE__ << "____iterator_____" << std::endl;
  std::cout << "\tnumber of volumes is " << vols.size() << std::endl;
  Range::iterator iter = vols.begin();
  
  for (unsigned i = 0; i<vols.size(); ++i, ++iter)
  {
      code = dagmc.measure_volume(*iter, volume_measure);
      std::cout << "\tvolume of entity " << i << " is " << volume_measure << std::endl;
  }

  // Julie's way - 
  std::cout << __FILE__ << ", " << __func__ << ":" << __LINE__ << "___entity_by_id___" << std::endl;
  std::cout << "\tnum_vols is " << num_vols << std::endl;
  EntityHandle entity = NULL;
  for (unsigned i = 0; i<num_vols; i++)
  {
      entity = DAG->entity_by_id(3, i);
      code = DAG->measure_volume(entity, volume_measure);
      std::cout << "\tvolume of entity " << i << " is " << volume_measure << std::endl;
  }
}


