########################################################################
#
# Julie Zachman
# December 14, 2012
# Rudimentary makefile to test compiling and linking of heterogenous sources 
# 
# Modified A Davis
# Jan 18th 2013

CXX = g++

# Modify MOAB to point to the local MOAB installation
MOAB = /data/opt/dagmc/moab
FLUDAG = /data/opt/fludag/DAGMC/FluDAG

# FLUPRO is an environment variable that should be set to the directory of the FLUKA installation
FLUKA_LIBS_LINK = -L${FLUPRO} -lflukahp

# MOAB_INCLUDES uses the ${MOAB} var to simplify compilation on different system
MOAB_INCLUDES  = -I${MOAB}/include
MOAB_LIBDIR    =   ${MOAB}/lib

# MOAB_FLAGS determine what libraries we wish to include. Note: netcdf does not seem to be needed
MOAB_FLAGS = -lMOAB -ldagmc
MOAB_LIBS_LINK = -L${MOAB_LIBDIR} -lm ${MOAB_FLAGS} 

CXXFLAGS = ${MOAB_FLAGS} 
FLUFOR = gfortran
OBJS = fluka_funcs.o userini.o


fludag: $(OBJS)
	${FLUPRO}/flutil/lfluka -o fludag -a MOAB -a dagmc -m userini -m fluka_funcs -L${FLUPRO} -L${MOAB_LIBDIR} 

fluka_funcs.o: 	fluka_funcs.cpp fluka_funcs.h
	$(CXX) $(CXXFLAGS) $(MOAB_INCLUDES) -c fluka_funcs.cpp 

userini.o: userini.f
	$(FLUPRO)/flutil/fff userini.f


clean:
	rm *.o
