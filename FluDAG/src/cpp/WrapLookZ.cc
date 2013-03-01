// FluDAG tag 
// WrapLookZ.hh - Sara Vanini 
//
// Wrapper for localisation of starting point of particle.
//
// modified 20/III/00: history initialization moved to ISVHWR
//////////////////////////////////////////////////////////////////

#include "DagWrappers.hh"
#include "DagWrapUtils.hh"

using namespace moab;

void lkwr(double& pSx, double& pSy, double& pSz,
          double* pV, const int& oldReg, const int& oldLttc,
          int& newReg, int& flagErr, int& newLttc)
{
  std::cerr << "======= LKWR =======" << std::endl;
  std::cerr << "oldReg is " << oldReg << std::endl;
  std::cerr << "position is " << pSx << " " << pSy << " " << pSz << std::endl; 


  const double xyz[] = {pSx, pSy, pSz}; // location of the particle (xyz)
  int is_inside = 0; // logical inside or outside of volume
  int num_vols = DAG->num_entities(3); // number of volumes

  for (int i = 1 ; i <= num_vols ; i++) // loop over all volumes
    {
      EntityHandle volume = DAG->entity_by_index(3, i); // get the volume by index
      // No ray history or ray direction.
      ErrorCode code = DAG->point_in_volume(volume, xyz, is_inside);

      // check for non error
      if(MB_SUCCESS != code) 
	{
	  std::cerr << "Error return from point_in_volume!" << std::endl;
	  flagErr = 1;
	  return;
	}
      
      if (is_inside == 1 )  // we are inside the cell tested
	{
	  newReg = i;
	  flagErr = i;
          //BIZZARLY - WHEN WE ARE INSIDE A VOLUME, BOTH, newReg has to equal flagErr
	  std::cerr << "newReg is " << newReg << std::endl;
	  return;
	}

    }

  std::cerr << "point is not in any volume" << std::endl;
  return;
}

