          ########################################################################
#
# Julie Zachman
# December 14, 2012
# Rudimentary makefile to test compiling and linking of heterogenous sources 
# 

# Modify MOAB to point to the local MOAB installation
# MOAB = /data/opt/dagmc/moab
MOAB = ${HOME}/dagmc_bld/MOAB
FLUDAG = ${HOME}/DAGMC/FluDAG
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


FSRCS = 

CCSRCS = 

OBJS = $(FSRCS:.f=.o) $(CCSRCS:.c=.o)
PRGS = $(FSRCS:.f=) $(CCSRCS:.c=)

.SUFFIXES:
.SUFFIXES: .f .c .o

.f.o:
	cd obj; make $<; cd ..

.c.o: 
	g++ -c $<

CPPFLAGS = ${MOAB_INCLUDES} ${FLUDAG_INCLUDES}
CXXFLAGS = -g

default: all

all: fluka_funcs.o userini.o

fludagW.o : cpp/fludagW.cpp
	g++ -c $(CPPFLAGS) $<
 
# Create a fortrain main linking object files compiled from fortran, C and C++

main_f.o: ../for/main_f.f
	gfortran -c $<

fluka_funcs.o: ../cpp/fluka_funcs.cpp
	g++ -c $(CPPFLAGS) $< 
	mv ../cpp/fluka_funcs.o .

userini.o: ../for/userini.f
	${FLUPRO}/flutil/fff ../for/userini.f
	mv ../for/userini.o .
clean:
	rm *.o

print-%:
	@echo $* = $($*)
