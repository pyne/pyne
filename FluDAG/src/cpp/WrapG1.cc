

///////////////////////////////////////////////////////////////////
//
// WrapG1.hh - Sara Vanini
//
// Wrapper for geometry tracking: returns approved step of 
// particle and alla variables that fluka G1 computes.
//
//
/////////////////////////////////////////////////////////////////////

#include "DagWrappers.hh"
#include "DagWrapUtils.hh"

using namespace moab;

void g1wr(double& pSx, double& pSy, double& pSz, double* pV,
          int& oldReg, const int& oldLttc, double& propStep,
          int& nascFlag, double& retStep, int& newReg,
          double& saf, int& newLttc, int& LttcFlag,
          double* sLt, int* jrLt)
{
  std::cerr<<"============= G1WR =============="<<std::endl;    

  std::cerr << "Position " << pSx << " " << pSy << " " << pSz << std::endl;
  std::cerr << "Direction vector " << pV[0] << " " << pV[1] << " " << pV[2] << std::endl;
  std::cerr << "Oldreg = " << oldReg << std::endl;
  std::cerr << "PropStep = " << propStep << std::endl;
  
  MBEntityHandle vol = DAG->entity_by_index(3,oldReg);
  //  MBEntityHandle prev = DAG->entity_by_index( 2, *jsu );
  double next_surf_dist;
  MBEntityHandle next_surf = 0;
  MBEntityHandle newvol = 0;

  double point[3] = {pSx,pSy,pSz};
  double dir[3]   = {pV[0],pV[1],pV[2]};  

  MBErrorCode result = DAG->ray_fire(vol, point, dir,next_surf, next_surf_dist );

  retStep = next_surf_dist;

  MBErrorCode rval = DAG->next_vol(next_surf,vol,newvol);
  
  newReg = DAG->index_by_handle(newvol);

  std::cerr << "newReg = " << newReg << std::endl;

  return;
}




