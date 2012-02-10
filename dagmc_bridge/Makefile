#MOAB_BASE_DIR = /local.hd/cnergg/sjackson/MOAB/cubit
MOAB_BASE_DIR = /path/to/MOAB/install

include ${MOAB_BASE_DIR}/lib/moab.make
#include ${MOAB_BASE_DIR}/lib/iMesh-Defs.inc

CXXFLAGS = -g -Wall -fPIC
#CXXFLAGS = -Wall
SONAME=dagmc_bridge.so
LDFLAGS= -shared -Wl,-soname,${SONAME} ${IMESH_LIBS} -L${MOAB_BASE_DIR}/lib -ldagmc -lMOAB
LDFLAGS+=-Wl,-rpath=${MOAB_BASE_DIR}/lib
CXXSOURCES = dagmc_bridge.cpp
CXXOBJS = ${CXXSOURCES:.cpp=.o}

${SONAME}: ${CXXOBJS} Makefile
	${CXX} ${CXXFLAGS} -o $@ ${CXXOBJS} ${LDFLAGS}

dagmc_bridge.o: dagmc_bridge.cpp dagmc_bridge.h

.cpp.o:
	${CXX} ${CXXFLAGS} ${MOAB_INCLUDES} -o $@ -c $<

clean:
	rm -rf ${SONAME} *.o
