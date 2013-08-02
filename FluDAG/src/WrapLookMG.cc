
// FluDag tag 

///////////////////////////////////////////////////////////////////
//
// WrapLookMG.hh - Sara Vanini 26/X/99
//
// Wrapper for localisation of particle for magnetic field tracking 
//
// modified 13/IV/00: check if point belongs to one of the lattice 
// histories stored in jrLtGeant 
//
///////////////////////////////////////////////////////////////////


#include "DagWrappers.hh"
#include "DagWrapUtils.hh"

void lkmgwr(double& pSx, double& pSy, double& pSz,
            double* pV, const int& oldReg, const int& oldLttc,
	    int& flagErr, int& newReg, int& newLttc)
{
    std::cerr<<"============= LKMGWR =============="<<std::endl;
    return;
}

