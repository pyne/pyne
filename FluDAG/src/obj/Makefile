########################################################################
#
# Julie Zachman
# March 25, 2013
# Compile object files
# 

# Makefile variables MOAB and FLUDAG should be modified in makefile.local
# to point to the local FluDAG and MOAB installation
include ../makefile.local
FLUDAG_SOURCE = ${FLUDAG}/src

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
CXXFLAGS += -m64 -pedantic

TARGETS = obj 

wrapobjects = WrapInit.o \
	WrapLookDB.o \
	WrapSavHist.o \
	WrapLookZ.o \
	WrapDN.o \
	WrapReg2Name.o \
	WrapLookMG.o \
	WrapMag.o \
	WrapLookFX.o \
	WrapReg.o \
	WrapIniHist.o \
	WrapIncrHist.o \
	WrapFlgfwr.o \
	WrapG1RT.o

default: obj

all: obj 

# Look in directory arg2 for file-pattern arg1
vpath %.cpp ../cpp
vpath %.cc  ../cpp
vpath %.f   ../for
vpath %.c   ../c
vpath %.o   . 
vpath %.h   ../include 
vpath %.hh  ../include


# Removed linking of createRegionAndmatFiles.o, FlukaMaterial.o
# mainfludag: mainFluDAG.o fluka_funcs.o $(wrapobjects) 
# gfortran -o $@ $^ $(MOAB_LIBS_LINK) $(FLUKA_LIBS_LINK) -lstdc++
# mv mainfludag ..


%.o: %.cpp
	g++ -c $(CXXFLAGS) $< 

obj: $(wrapobjects) mainFluDAG.o fluka_funcs.o
