########################################################################
#
# Julie Zachman
# December 14, 2012
# Rudimentary makefile to test compiling and linking of heterogenous sources 
# 

# Modify MOAB to point to the local MOAB installation
# MOAB = /data/opt/dagmc/moab
MOAB = ${HOME}/dagmc_bld/MOAB
FLUDAG = ${HOME}/storage/DAGMC/FluDAG
FLUDAG_SOURCE = ${FLUDAG}/source

# FLUPRO is an environment variable that should be set to the directory of the FLUKA installation
FLUKA_LIBS_LINK = -L${FLUPRO} -lflukahp

# Local includes
FLUDAG_INCLUDES = -I${FLUDAG_SOURCE}/include

# MOAB_INCLUDES uses the ${MOAB} var to simplify compilation on different systems
MOAB_INCLUDES  = -I${MOAB}/include
MOAB_LIBDIR    =   ${MOAB}/lib

# MOAB_FLAGS determine what libraries we wish to include. Note: netcdf does not seem to be needed
MOAB_FLAGS = -lMOAB -ldagmc
MOAB_LIBS_LINK = -L${MOAB_LIBDIR} -lm ${MOAB_FLAGS} 

CXXFLAGS = ${MOAB_INCLUDES} ${FLUDAG_INCLUDES} # Modify this to the name of the directory from which myfluka is run
RUNDIR = ./userini_example

TARGETS = mainfludag mainf

default: all

all: myfluka mainfludag

# Look in directory arg2 for file-pattern arg1
vpath %.cpp cpp
vpath %.f   for
vpath %.c   c
vpath %.o   obj

mainfludag: mainFluDAG.o fluka_funcs.o readVol.o createRegionAndMatFiles.o
	gfortran -o $@ $^ $(MOAB_LIBS_LINK) $(FLUKA_LIBS_LINK) -lstdc++
# Try this as well - 
#	 g++ -o $@ $^ $(MOAB_LIBS_LINK) $(FLUKA_LIBS_LINK) 
	
myfludag: 
# Create a fortran main linking object files compiled from fortran, C and C++
mainf: main_f.o fluka_funcs.o
#	gfortran -o $@ obj/main_f.o obj/fluka_funcs.o $(MOAB_LIBS_LINK) $(FLUKA_LIBS_LINK) -lstdc++
	@echo $^
	gfortran -o $@ $^ $(MOAB_LIBS_LINK) $(FLUKA_LIBS_LINK) -lstdc++

# How to compile and link a user routine: FM.pdf page 341.  See also "$FLUPRO/flutil/lfluka -h"
# This compiles a fluka main that will run the code in userini.f.  Any calls userini.f makes to 
# other libraries is linked in at the end.  In this case, userini makes a call into fluka_funcs.
# Important:  userini must be so named.
myfluka: userini.o fluka_funcs.o
	${FLUPRO}/flutil/lfluka -o $@ -a MOAB -a dagmc -m fluka obj/userini.o -L ${FLUPRO} -L ${MOAB_LIBDIR} obj/fluka_funcs.o

main_f.o : main_f.f
	gfortran -c $<
	mv $@ obj

userini.o : userini.f
	${FLUPRO}/flutil/fff $<
	mv for/$@ obj
	 
fluka_funcs.o : fluka_funcs.cpp
	g++ -c $(CXXFLAGS) $< 
	mv $@ obj

mainFluDAG.o : mainFluDAG.cpp 
	g++ -c $(CXXFLAGS) $< 
	mv $@ obj

readVol.o : readVol.cpp 
	g++ -c $(CXXFLAGS) $< 
	mv $@ obj

createRegionAndMatFiles.o : createRegionAndMatFiles.cpp
	g++ -c $(CXXFLAGS) $< 
	mv $@ obj


test.stl : test.h5m
	mbconvert $< $@

install : 
	cp ${TARGETS} ../bin

.PHONY : clean
clean:
	rm obj/*.o
	rm myfluka
print-%:
	@echo $* = $($*)
