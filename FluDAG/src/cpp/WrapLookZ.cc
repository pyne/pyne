// FluDAG tag 

///////////////////////////////////////////////////////////////////
//
// WrapLookZ.hh - Sara Vanini 
//
// Wrapper for localisation of starting point of particle.
//
// modified 20/III/00: history initialization moved to ISVHWR
//////////////////////////////////////////////////////////////////


#include "DagWrappers.hh"
#include "DagWrapUtils.hh"

using namespace moab;

#define VOL 3

void lkwr(double& pSx, double& pSy, double& pSz,
          double* pV, const int& oldReg, const int& oldLttc,
          int& newReg, int& flagErr, int& newLttc)
{
  std::cout << "======= LKWR =======" << std::endl;
  std::cout << "oldReg is " << oldReg << std::endl;

  const double xyz[] = {pSx, pSy, pSz};
  int is_inside;
  int num_vols = DAG->num_entities(3);
   
  for (int i = 1 ; i <= num_vols ; i++)
  {
    EntityHandle volume = DAG->entity_by_index(VOL, i);
    // No ray history or ray direction.
    ErrorCode code = DAG->point_in_volume(volume, xyz, is_inside, NULL, NULL);
    if(MB_SUCCESS != code) 
    {
        std::cout << "Error return from point_in_volume!" << std::endl;
        return;
    }
    if (is_inside) 
    {
	// Which is it?
        // newReg = DAG->id_by_entity(i);
        newReg = i;
        std::cout << "newReg is " << newReg << std::endl;
        return;
    }
  }
  std::cout << "newReg is not found" << std::endl;
  
}

/*

dagmcchkcel_(double *uuu,double *vvv,double *www,double *xxx,
                  double *yyy,double *zzz, int *i1, int *j)
{


#ifdef TRACE_DAGMC_CALLS
  std::cout<< " " << std::endl;
  std::cout<< "chkcel: vol=" << DAG->id_by_index(3,*i1) << " xyz=" << *xxx
           << " " << *yyy << " " << *zzz << std::endl;
  std::cout<< "      : uvw = " << *uuu << " " << *vvv << " " << *www << std::endl;
#endif

  int inside;
  MBEntityHandle vol = DAG->entity_by_index( 3, *i1 );
  double xyz[3] = {*xxx, *yyy, *zzz};
  double uvw[3] = {*uuu, *vvv, *www};
  MBErrorCode rval = DAG->point_in_volume( vol, xyz, inside, uvw );

  if (MB_SUCCESS != rval) {
    std::cerr << "DAGMC: failed in point_in_volume" <<  std::endl;
    exit(EXIT_FAILURE);
  }

  if (MB_SUCCESS != rval) *j = -2;
  else
    switch (inside)
      {
      case 1:
        *j = 0; // inside==  1 -> inside volume -> j=0
        break;
      case 0:
        *j = 1; // inside== 0  -> outside volume -> j=1
        break;
      case -1:
        *j = 1; // inside== -1 -> on boundary -> j=1 (assume leaving volume)
        break;
      default:
        std::cerr << "Impossible result in dagmcchkcel" << std::endl;
        exit(EXIT_FAILURE);
      }

#ifdef TRACE_DAGMC_CALLS
  std::cout<< "chkcel: j=" << *j << std::endl;
#endif

}
*/
